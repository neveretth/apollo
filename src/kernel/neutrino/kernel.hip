#define __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>

#include "../../defs.h"

extern "C" __device__ real_t compute_next_timestep_neutrino(real_t* E_O,
                                     real_t* E_D,
                                     real_t dt_old, int i, int num_groups) {
    real_t dt_new = (real_t)0.0;
    // if (i == 0) {
        real_t E_R = (real_t)0.0;
        for (uint k = 0; k < num_groups; k++) {
            // E_R += E_O[k] / E_D[k];
            // Atomic XOR, largest log2 OOM?
            E_R = max(E_R, E_O[k] / E_D[k]);
        }

        if (E_R < dt_tol_low) {
            dt_new = dt_old * 0.9 * pow(E_R, (real_t)-1.0);
        } else if (E_R < dt_tol_high) {
            dt_new = dt_old * 0.9 * pow(E_R, (real_t)-0.5);
        } else {
            dt_new = dt_old;
        }

        dt_new = max(min(dt_new, (real_t)2.0 * dt_old),
                     (real_t)0.5 * dt_old);
    // }
    return dt_new;
}

extern "C" __device__ void QSS1(real_t* N_p, real_t* Alpha,
          real_t* F, real_t* k, real_t dt,
          real_t* Nold, int i) {
    N_p[i] = Nold[i];
    real_t r;

    r = 1.0 / (k[i] * dt);

    Alpha[i] = ((160 * pow(r, (real_t)3)) +
                (60 * pow(r, (real_t)2)) + (11 * r) + 1) /
               ((360 * pow(r, (real_t)3)) +
                (60 * pow(r, (real_t)2)) + (12 * r) + 1);

    N_p[i] += dt * (F[i] - (k[i] * Nold[i])) /
              (1.0 + (Alpha[i] * k[i] * dt));
}

extern "C" __device__ void QSS2(real_t* Nnew, real_t* F0,
          real_t* Fp, real_t* k0,
          real_t* kp, real_t* Alpha0,
          real_t dt, real_t* Nold, int i) {
    Nnew[i] = Nold[i];

    real_t Ft;
    real_t kBAR;
    real_t rBAR;
    real_t AlphaBAR;

    kBAR = 0.5 * (kp[i] + k0[i]);
    rBAR = 1.0 / (kBAR * dt);
    AlphaBAR =
        ((160 * pow(rBAR, (real_t)3)) +
         (60 * pow(rBAR, (real_t)2)) + (11 * rBAR) + 1) /
        ((360 * pow(rBAR, (real_t)3)) +
         (60 * pow(rBAR, (real_t)2)) + (12 * rBAR) + 1);
    Ft = AlphaBAR * Fp[i] + (1.0 - AlphaBAR) * F0[i];
    Nnew[i] += dt * (Fp[i] - kBAR * Nold[i]) /
               (1.0 + (Alpha0[i] * kBAR * dt));
}

extern "C" __device__ void compute_rates(real_t* F, real_t* k,
                   real_t* R_In, real_t* R_Out, real_t* N,
                   real_t* dV, int i, int i1d, int num_groups) {
    k[i] = (real_t)0.0;
    F[i] = (real_t)0.0;
    for (int j = 0; j < num_groups; j++) {
        real_t temp = R_In[i1d + j] * N[j] * dV[j];
        F[i] += temp;
        k[i] += temp + R_Out[i1d + j] * (1.0 - N[j]) * dV[j];
    }
}

extern "C" __global__ void qss_hip_neunet_kernel(real_t* real_t_val, int* int_val,
                               real_t* R_In, real_t* R_Out, real_t* eC, real_t* dV,
                               real_t* N_Eq, real_t* Nold) {
    const int i = threadIdx.x;

    int num_groups;
    real_t t;
    real_t dt;
    real_t G_A;
    real_t G_B;
    real_t G_C;
    real_t tol;
    real_t EpsA;
    real_t EpsR;
    real_t tolC;
    real_t ERR;
    real_t t_end;
    // if (i == 0) {
        // t = real_t_val[0];
        t = 0; // This should just be zero.
        dt = real_t_val[1];
        G_A = real_t_val[2];
        G_B = real_t_val[3];
        G_C = real_t_val[4];
        tol = real_t_val[5];
        EpsA = real_t_val[6];
        EpsR = real_t_val[6];
        tolC = real_t_val[8];
        ERR = real_t_val[9];
        t_end = real_t_val[10];
        num_groups = int_val[0];
    // }

    const int i1d = num_groups * i;

    __syncthreads();

    bool done = false;
    bool reStep = false;
    int cycle = 0;
    int true_cycle = 0;
    real_t dt_new = dt;

    real_t* N_0 = new real_t[num_groups];
    real_t* kp = new real_t[num_groups];
    real_t* Fp = new real_t[num_groups];
    real_t* Np = new real_t[num_groups];
    real_t* F0 = new real_t[num_groups];
    real_t* k0 = new real_t[num_groups];
    real_t* F1 = new real_t[num_groups];
    real_t* k1 = new real_t[num_groups];
    real_t* N_1 = new real_t[num_groups];
    real_t* N_2 = new real_t[num_groups];
    real_t* Alpha0 = new real_t[num_groups];
    real_t* Error_Observed = new real_t[num_groups];
    real_t* Error_Desired = new real_t[num_groups];
    real_t* TEMP = new real_t[num_groups];
    real_t* tempN1 = new real_t[num_groups];
    real_t* Nnew = new real_t[num_groups];
    real_t* CN = new real_t[num_groups];
    
    __syncthreads();

    N_0[i] = G_A * exp(-0.5 * pow((eC[i] - G_B) / G_C, 2));

    Nold[i] = N_0[i];

    __syncthreads();

    while (!done) {
        true_cycle++;

        if (!reStep) {
            cycle++;
        }

        // It goes without saying that this is broken.
        // But the point is to get a concept of computational expense...
        // This is approximately right to make the time accurate.
        if (true_cycle > 300) {
            return;
        }

        {
            Nnew[i] = 0;
          
            dt = dt_new;

            compute_rates(F0, k0, R_In, R_Out, Nold, dV, i,
                          i1d, num_groups);

            QSS1(Np, Alpha0, F0, k0, 0.5 * dt, Nold, i);
            __syncthreads();

            compute_rates(Fp, kp, R_In, R_Out, Np, dV, i,
                          i1d, num_groups);

            QSS2(N_1, F0, Fp, k0, kp, Alpha0, 0.5 * dt, Nold, i);
            __syncthreads();

            compute_rates(F1, k1, R_In, R_Out, Nold, dV, i,
                          i1d, num_groups);

            QSS1(Np, Alpha0, F1, k1, 0.5 * dt, N_1, i);
            __syncthreads();

            compute_rates(Fp, kp, R_In, R_Out, Np, dV, i,
                          i1d, num_groups);

            QSS2(N_1, F0, Fp, k0, kp, Alpha0, 0.5 * dt, N_1, i);
            __syncthreads();

            QSS1(Np, Alpha0, F0, k0, dt, Nold, i);
            __syncthreads();

            compute_rates(Fp, kp, R_In, R_Out, Np, dV, i,
                          i1d, num_groups);

            QSS2(N_2, F0, Fp, k0, kp, Alpha0, dt, Nold, i);
            __syncthreads();
            
            Nnew[i] = N_1[i];
            __syncthreads();

            Error_Observed[i] = fabs(N_1[i] - N_2[i]);
            Error_Desired[i] = EpsA + EpsR * fabs(N_1[i]);
            __syncthreads();
            
            dt_new = compute_next_timestep_neutrino(
                Error_Observed, Error_Desired, dt, i, num_groups);
        }

        Nold[i] = Nnew[i];
        t += dt;

        if ((t - t_end) >= 0) {
            done = true;
            break;
        }
    }

    if (i == 0) {
        real_t_val[0] = t;
        real_t_val[1] = dt;
    }
    
    Nold[i] = Nold[i];

    delete [] N_0;
    delete [] kp;
    delete [] Fp;
    delete [] Np;
    delete [] F0;
    delete [] k0;
    delete [] F1;
    delete [] k1;
    delete [] N_1;
    delete [] N_2;
    delete [] Alpha0;
    delete [] Error_Observed;
    delete [] Error_Desired;
    delete [] TEMP;
    delete [] tempN1;
    delete [] Nnew;
    delete [] CN;
    
    return;
}
