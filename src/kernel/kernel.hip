// Omit:
// #include "kernel.h"
// But access to these functions is still granted through kernel.h ...
// This is profoundly stupid, but HIP has forced my hand.
// Let me do C! Some people use C!
#include <hip/hip_runtime.h>

#define THIRD 0.333333333333333

extern "C" __global__ void vector_mult_kernel(float* A, float* B, float* C) {
    int tid = threadIdx.x;
    C[tid] = A[tid] * B[tid];
}

extern "C" __device__ int check_asy(float Fminus, float Y, float dt) {
    if (Y > 0.0f && Fminus * dt / Y > 1.0f) {
        return 1;
    } else {
        return 0;
    }
}

extern "C" __device__ float asymptotic_update(float Fplus, float Fminus, float Y, float dt) {
    return (Y + Fplus * dt) / (1.0f + Fminus * dt / Y); // Sophia He formula
}

extern "C" __device__ float euler_update(float FplusSum, float FminusSum, float dt) {
    return (FplusSum - FminusSum) * dt;
}

extern "C" __device__ float compute_timestep(float prevdt, float t, float tmax) {
    float dt;
    if (t == 0.0f) {
        dt = 1.0e-20;
    } else {
        dt = 0.1f * t;
    }

    // Prevent final integration step from overstepping tmax
    // if(t+dt > tmax) dt = tmax - t;

    return dt;
}

extern "C" __device__ float compute_keff(float Fminus, float Y) {
    if (Y > 0) {
        return Fminus / Y;
    } else {
        return 0.0f;
    }
}

extern "C" __global__ void integration_kernel_hip(float* P0, float* P1, float* P2,
        float* P3, float* P4, float* P5, float* P6, float* Prefac, float* Q,
        float* Rate, float* Flux, float* Fplus, float* Fminus, float* FplusFac,
        float* FminusFac, float* FplusSum, float* FminusSum, int* FplusMax,
        int* FminusMax, int* MapFplus, int* MapFminus, float* Y,
        int* NumReactingSpecies, int* Reactant1, int* Reactant2,
        int* Reactant3, int* int_val, float* float_val) {

    int integration_steps = 0;
    
    int number_species = int_val[0];
    int number_reactions = int_val[1];
    int f_plus_total = int_val[2];
    int f_minus_total = int_val[3];
    int halt = int_val[4];
    float t9 = float_val[0];
    float t_max = float_val[1];
    float dt_init = float_val[2];

    int i = threadIdx.x;

    // Compute the temperature-dependent factors for the rates.
    float t93 = powf(t9, THIRD);
    float t1 = 1 / t9;
    float t2 = 1 / t93;
    float t3 = t93;
    float t4 = t9;
    float t5 = t93 * t93 * t93 * t93 * t93;
    float t6 = logf(t9);

    Rate[i] =
        Prefac[i] * expf(P0[i] + t1 * P1[i] + t2 * P2[i] + t3 * P3[i] +
                t4 * P4[i] + t5 * P5[i] + t6 * P6[i]);

    __syncthreads();

    /*
     * Begin the time integration from t=0 to tmax. Rather than t=0 we start at
     * some very small value of t.
     */
    float t = 1.0e-16;      // The current integration time
    float dt = dt_init;     // The current integration timestep
    float prevdt = dt_init; // The integration timestep from the previous step

    // Main time integration loop
    while (t < t_max) {
        // Compute the fluxes from the previously-computed rates and the current
        // abundances
        if (i < number_reactions) {
            int nr = NumReactingSpecies[i];
            switch (nr) {
                case 1:
                    Flux[i] = Rate[i] * Y[*(Reactant1 + i)];
                    break;
                case 2:
                    Flux[i] = Rate[i] * Y[*(Reactant1 + i)] * Y[*(Reactant2 + i)];
                    break;
                case 3:
                    Flux[i] = Rate[i] * Y[*(Reactant1 + i)] * Y[*(Reactant2 + i)] *
                        Y[*(Reactant3 + i)];
                    break;
            }
        }
        __syncthreads();
    
        // Populate the F+ and F- arrays in parallel from the master Flux array
        if (i < f_plus_total) {
            Fplus[i] = FplusFac[i] * Flux[MapFplus[i]];
        }
        __syncthreads();
        
        if (i < f_minus_total) {
            Fminus[i] = FminusFac[i] * Flux[MapFminus[i]];
        }
        __syncthreads();

        if (i < number_species) {
            // Partially serial Sum F+
            int minny = 0;
            if (i > 0) {
                minny = FplusMax[i - 1] + 1;
            }
            FplusSum[i] = 0.0f;
            for (int j = minny; j <= FplusMax[i]; j++) {
                FplusSum[i] += Fplus[j];
            }

            // Partially serial Sum F-
            minny = 0;
            if (i > 0) {
                minny = FminusMax[i - 1] + 1;
            }
            FminusSum[i] = 0.0f;
            for (int j = minny; j <= FminusMax[i]; j++) {
                FminusSum[i] += Fminus[j];
            }
        }
        __syncthreads();

        /*
           Now use the fluxes to update the populations in parallel for this
           timestep For now we shall assume the asymptotic method. We determine
           whether each isotope satisfies the asymptotic condition. If it does we
           update with the asymptotic formula. If not, we update numerically using
           the forward Euler formula.
         */
        if (i < number_species) {
            if (check_asy(Fminus[i], Y[i], dt) == 1) {
                Y[i] = asymptotic_update(FplusSum[i], FminusSum[i], Y[i], dt);
            } else {
                Y[i] += euler_update(FplusSum[i], FminusSum[i], dt);
            }
        }
        __syncthreads();

        // Increment the integration time and set the new timestep
        t += dt;
        integration_steps++;
        prevdt = dt;
        dt = compute_timestep(prevdt, t, t_max);
        __syncthreads();

        // Temporary diagnostic halt
        if (integration_steps >= halt) {
            break;
        }
    }
}
