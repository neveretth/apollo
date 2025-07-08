#include "kernel.h"

#include <math.h>
#include <stdlib.h>

#define THIRD 0.3333333333333333

int tnn_integration_kernel(float* P0, float* P1, float* P2, float* P3,
                           float* P4, float* P5, float* P6, float* Prefac,
                           float* Q, float* Rate, float* Flux, float* Fplus,
                           float* Fminus, float* FplusFac, float* FminusFac,
                           float* FplusSum, float* FminusSum, int* FplusMax,
                           int* FminusMax, int* MapFplus, int* MapFminus,
                           float* Y, int* NumReactingSpecies, int* Reactant1,
                           int* Reactant2, int* Reactant3, int number_species,
                           int number_reactions, int f_plus_total,
                           int f_minus_total, float t9, float t_max,
                           float dt_init, int halt) {

    int integration_steps = 0;

    // Compute the temperature-dependent factors for the rates.
    float t93 = powf(t9, THIRD);
    float t1 = 1 / t9;
    float t2 = 1 / t93;
    float t3 = t93;
    float t4 = t9;
    float t5 = t93 * t93 * t93 * t93 * t93;
    float t6 = logf(t9);

    for (int i = 0; i < number_reactions; i++) {
        Rate[i] =
            Prefac[i] * expf(P0[i] + t1 * P1[i] + t2 * P2[i] + t3 * P3[i] +
                             t4 * P4[i] + t5 * P5[i] + t6 * P6[i]);
    }

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
        for (int i = 0; i < number_reactions; i++) {
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

        // Populate the F+ and F- arrays in parallel from the master Flux array
        for (int i = 0; i < f_plus_total; i++) {
            Fplus[i] = FplusFac[i] * Flux[MapFplus[i]];
        }

        for (int i = 0; i < f_minus_total; i++) {
            Fminus[i] = FminusFac[i] * Flux[MapFminus[i]];
        }

        for (int i = 0; i < number_species; i++) {
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

        /*
        Now use the fluxes to update the populations in parallel for this
        timestep For now we shall assume the asymptotic method. We determine
        whether each isotope satisfies the asymptotic condition. If it does we
        update with the asymptotic formula. If not, we update numerically using
        the forward Euler formula.
        */
        for (int i = 0; i < number_species; i++) {
            if (check_asy(Fminus[i], Y[i], dt) == 1) {
                Y[i] = asymptotic_update(FplusSum[i], FminusSum[i], Y[i], dt);
            } else {
                Y[i] += euler_update(FplusSum[i], FminusSum[i], dt);
            }
        }

        // Increment the integration time and set the new timestep
        t += dt;
        integration_steps++;
        prevdt = dt;
        dt = compute_timestep(prevdt, t, t_max);

        // Temporary diagnostic halt
        if (integration_steps >= halt) {
            break;
        }
    }

    return EXIT_SUCCESS;
}

int check_asy(float Fminus, float Y, float dt) {
    if (Y > 0.0f && Fminus * dt / Y > 1.0f) {
        return 1;
    } else {
        return 0;
    }
}

float asymptotic_update(float Fplus, float Fminus, float Y, float dt) {
    return (Y + Fplus * dt) / (1.0f + Fminus * dt / Y); // Sophia He formula
}

float euler_update(float FplusSum, float FminusSum, float dt) {
    return (FplusSum - FminusSum) * dt;
}

float compute_timestep(float prevdt, float t, float tmax) {
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

float compute_keff(float Fminus, float Y) {
    if (Y > 0) {
        return Fminus / Y;
    } else {
        return 0.0f;
    }
}

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
    float dt_ea;
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
    float* exp_kdt = malloc(sizeof(float) * n_g);
    float* inv_k = malloc(sizeof(float) * n_g);
    float* alpha_0 = malloc(sizeof(float) * n_g);
    float* error_desired = malloc(sizeof(float) * n_g);
    float* error_observed = malloc(sizeof(float) * n_g);

    float* f0_k0;

    while (done == 0) {
        true_cycle++;

        if (!restep) {
            cycle++;
        }

        if (true_cycle >= halt) {
            return EXIT_SUCCESS;
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
            done = true;
            break;
        }
    }

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
        float Ft;
        float kBAR;
        float rBAR;
        float AlphaBAR;

        kBAR = 0.5 * (kp[i] + k0[i]);
        rBAR = 1.0 / (kBAR * dt);
        AlphaBAR =
            ((160 * pow(rBAR, 3)) + (60 * pow(rBAR, 2)) + (11 * rBAR) + 1) /
            ((360 * pow(rBAR, 3)) + (60 * pow(rBAR, 2)) + (12 * rBAR) + 1);
        Ft = AlphaBAR * fp[i] + (1.0 - AlphaBAR) * f0[i];
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
    for (size_t i = 0; i < n_g; i++) {
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

int hydro_integration_kernel(float* temp, float* density, float volume, float h,
                             float dt, float t_end, int dim) {
    float* delta_temp = malloc(dim * sizeof(float));
    float c = 1e8;  // Placeholder contribution value.
    float area = 1; // Placeholder area of interaction.
    float ntd = 0;  // Nuclear burning temp diff (assume negligible for now)

    float t = 0;

    // Does one integration (single step if t_end = 0)
    do {
        for (int i = 0; i < dim; i++) {
            float pdvh = 1 / (density[i] * volume * c);
            float temp_diff = -(2 * temp[i]);

            // Binary expression becomes 1 or 0. Faster than if statements.
            temp_diff += temp[i - (i > 0)];
            temp_diff += temp[i + (i < dim - 1)];

            delta_temp[i] = dt * (pdvh * (ntd + (h * area * temp_diff)));
        }
        // Yes, these need to be separate.
        for (int i = 0; i < dim; i++) {
            temp[i] += delta_temp[i];
        }
        t += dt;
    } while (t < t_end);

    free(delta_temp);

    return EXIT_SUCCESS;
}

int flat_hydro_integration_kernel(float** temp, float** density, float volume,
                                  float h, float dt, float t_end, int* dim) {
    float** delta_temp = malloc(dim[0] * sizeof(float*));
    for (int i = 0; i < dim[0]; i++) {
        delta_temp[i] = malloc(dim[1] * sizeof(float));
    }

    float c = 1e8;  // Placeholder contribution value.
    float area = 1; // Placeholder area of interaction.
    float ntd = 0;  // Nuclear burning temp diff (assume negligible for now)

    float t = 0;

    // Does one integration (single step if t_end = 0)
    do {
        for (int i = 0; i < dim[0]; i++) {
            for (int j = 0; j < dim[1]; j++) {
                float pdvh = 1 / (density[i][j] * volume * c);
                float temp_diff = -(4 * temp[i][j]);

                // As with the linear integration, towards idx=0 is positive
                temp_diff += temp[i - (i > 0)][j];
                temp_diff += temp[i + (i < dim[0] - 1)][j];
                temp_diff += temp[i][j - (j > 0)];
                temp_diff += temp[i][j + (j < dim[1] - 1)];

                delta_temp[i][j] = dt * (pdvh * (ntd + (h * area * temp_diff)));
            }
        }
        for (int i = 0; i < dim[0]; i++) {
            for (int j = 0; j < dim[1]; j++) {
                temp[i][j] += delta_temp[i][j];
            }
        }
        t += dt;
    } while (t < t_end);

    for (int i = 0; i < dim[0]; i++) {
        free(delta_temp[i]);
    }
    free(delta_temp);

    return EXIT_SUCCESS;
}
