#include "kernel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

int neunet_integration_kernel(float** rate_in, float** rate_out, float* n_old,
                              float* ec, float* dv, float dt, float t_end,
                              float EpsA, float EpsR, float g_a, float g_b,
                              float g_c, int n_g, int halt) {
    float* n_0 = malloc(sizeof(float) * n_g);

    for (int i = 0; i < n_g; i++) {
        n_0[i] = g_a * exp(-0.5 * pow((ec[i] - g_b) / g_c, 2));
    }

    for (int i = 0; i < n_g; i++) {
        n_old[i] = n_0[i];
    }

    int done = 0;
    int restep = 0;
    int cycle = 0;
    int true_cycle = 0;
    float dt_new = dt;
    float t = 0;

    float* n_new = malloc(sizeof(float) * n_g);
    float* n_1 = malloc(sizeof(float) * n_g);
    float* n_2 = malloc(sizeof(float) * n_g);
    float* k0 = malloc(sizeof(float) * n_g);
    float* k1 = malloc(sizeof(float) * n_g);
    float* f0 = malloc(sizeof(float) * n_g);
    float* f1 = malloc(sizeof(float) * n_g);
    float* fp = malloc(sizeof(float) * n_g);
    float* kp = malloc(sizeof(float) * n_g);
    float* np = malloc(sizeof(float) * n_g);
    float* alpha_0 = malloc(sizeof(float) * n_g);
    float* error_desired = malloc(sizeof(float) * n_g);
    float* error_observed = malloc(sizeof(float) * n_g);

    while (done == 0) {
        true_cycle++;

        if (!restep) {
            cycle++;
        }

        if (true_cycle >= halt) {
            goto exit;
        }

        memset(n_new, 0, n_g);

        dt = dt_new;

        compute_rates(f0, k0, rate_in, rate_out, n_old, dv, n_g);

        QSS1(np, alpha_0, f0, k0, 0.5 * dt, n_old, n_g);

        compute_rates(fp, kp, rate_in, rate_out, np, dv, n_g);

        QSS2(n_1, f0, fp, k0, kp, alpha_0, 0.5 * dt, n_old, n_g);

        compute_rates(f1, k1, rate_in, rate_out, n_old, dv, n_g);

        QSS1(np, alpha_0, f1, k1, 0.5 * dt, n_1, n_g);

        compute_rates(fp, kp, rate_in, rate_out, np, dv, n_g);

        QSS2(n_1, f0, fp, k0, kp, alpha_0, 0.5 * dt, n_1, n_g);

        QSS1(np, alpha_0, f0, k0, dt, n_old, n_g);

        compute_rates(fp, kp, rate_in, rate_out, np, dv, n_g);

        QSS2(n_2, f0, fp, k0, kp, alpha_0, dt, n_old, n_g);

        for (int i = 0; i < n_g; i++) {
            n_new[i] = n_1[i];
        }

        for (int i = 0; i < n_g; i++) {
            error_observed[i] = fabs(n_1[i] - n_2[i]);
            error_desired[i] = EpsA + EpsR * fabs(n_1[i]);
        }

        dt_new = compute_next_timestep(error_observed, error_desired, dt, n_g);

        for (int i = 0; i < n_g; i++) {
            n_old[i] = n_new[i];
        }

        t += dt;

        if ((t - t_end) >= 0) {
            done = 1;
            break;
        }
    }

exit:
    free(n_0);
    free(n_new);
    free(n_1);
    free(n_2);
    free(k0);
    free(k1);
    free(f0);
    free(f1);
    free(fp);
    free(kp);
    free(np);
    free(alpha_0);
    free(error_desired);
    free(error_observed);

    return EXIT_SUCCESS;
}

void QSS1(float* N_p, float* Alpha, float* F, float* k, float dt, float* n_old,
          int n_g) {
    for (int i = 0; i < n_g; i++) {
        N_p[i] = n_old[i];
        float r;

        r = 1.0 / (k[i] * dt);

        Alpha[i] =
            ((160 * pow(r, (float)3)) + (60 * pow(r, (float)2)) + (11 * r) +
             1) /
            ((360 * pow(r, (float)3)) + (60 * pow(r, (float)2)) + (12 * r) + 1);

        N_p[i] +=
            dt * (F[i] - (k[i] * n_old[i])) / (1.0 + (Alpha[i] * k[i] * dt));
    }
}

void QSS2(float* Nnew, float* f0, float* fp, float* k0, float* kp,
          float* alpha_0, float dt, float* n_old, int n_g) {
    for (int i = 0; i < n_g; i++) {
        Nnew[i] = n_old[i];
        // float Ft;
        float kBAR;
        float rBAR;
        float AlphaBAR;

        kBAR = 0.5 * (kp[i] + k0[i]);
        rBAR = 1.0 / (kBAR * dt);
        AlphaBAR =
            ((160 * pow(rBAR, 3)) + (60 * pow(rBAR, 2)) + (11 * rBAR) + 1) /
            ((360 * pow(rBAR, 3)) + (60 * pow(rBAR, 2)) + (12 * rBAR) + 1);
        // Ft = AlphaBAR * fp[i] + (1.0 - AlphaBAR) * f0[i];
        Nnew[i] +=
            dt * (fp[i] - kBAR * n_old[i]) / (1.0 + (alpha_0[i] * kBAR * dt));
    }
}

void compute_rates(float* F, float* k, float** R_In, float** R_Out, float* N,
                   float* dv, int n_g) {

    for (int i = 0; i < n_g; i++) {
        k[i] = 0.0;
        F[i] = 0.0;
        for (int j = 0; j < n_g; j++) {
            F[i] += R_In[i][j] * N[j] * dv[j];
            k[i] +=
                R_In[i][j] * N[j] * dv[j] + R_Out[i][j] * (1.0 - N[j]) * dv[j];
        }
    }
}

// "Temporary" for MAX/MIN
#include <sys/param.h>

float compute_next_timestep(const float* E_O, const float* E_D, float dt_old,
                            int n_g) {

    float E_R = 0.0;
    for (int i = 0; i < n_g; i++) {
        E_R = MAX(E_R, E_O[i] / E_D[i]);
    }

    float dt_new;
    if (E_R > 0.1) {
        dt_new = dt_old * 0.9 * pow(E_R, -1.0);
    } else if (E_R < 0.5) {
        dt_new = dt_old * 0.9 * pow(E_R, -0.5);
    } else {
        dt_new = dt_old;
    }

    dt_new = MAX(MIN(dt_new, 2.0 * dt_old), 0.5 * dt_old);

    return dt_new;
}

#define BOLTZMANN_CONSTANT 1.380658e-16

// This really highlights the weakness of storing the neutrino network data like
// this. A rework could be extremely beneficial for parallelization.
int neunet_data_preprocess(struct neunet**** neunet, struct rt_hydro_mesh* mesh,
                           struct simulation_properties sim_prop) {
    for (int i = 0; i < sim_prop.resolution[0]; i++) {
        for (int j = 0; j < sim_prop.resolution[1]; j++) {
            for (int k = 0; k < sim_prop.resolution[2]; k++) {
                neunet[i][j][k]->f->kt =
                    mesh->temp[i][j][k] * BOLTZMANN_CONSTANT;
                neunet[i][j][k]->f->rho = mesh->density[i][j][k];
                neunet[i][j][k]->f->t_end = mesh->dt;
            }
        }
    }
    return EXIT_SUCCESS;
}
