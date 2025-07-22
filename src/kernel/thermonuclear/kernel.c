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
