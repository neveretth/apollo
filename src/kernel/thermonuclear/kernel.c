#include "kernel.h"

#include <math.h>
#include <stdlib.h>

#define THIRD 0.3333333333333333
#define ECON 9.5768e17 // Convert MeV/nucleon/s to erg/g/s

// This code has a remarkably bad habit of nan/inf-ing out.
int tnn_integration_kernel(
    real_t* real_t_val, int* int_val, real_t* P0, real_t* P1, real_t* P2,
    real_t* P3, real_t* P4, real_t* P5, real_t* P6, real_t* Prefac, real_t* Q,
    real_t* Rate, real_t* Flux, real_t* Fplus, real_t* Fminus, real_t* FplusFac,
    real_t* FminusFac, real_t* FplusSum, real_t* FminusSum, int* FplusMax,
    int* FminusMax, int* MapFplus, int* MapFminus, real_t* Y,
    int* NumReactingSpecies, int* Reactant1, int* Reactant2, int* Reactant3,
    int number_species, int number_reactions, int f_plus_total,
    int f_minus_total, real_t t9, real_t t_max, real_t dt_init) {

    int integration_steps = 0;

    // Placeholder transfer for tmp vars, this will be implemented properly in
    // a future version.
    t9 = real_t_val[0];

    // Compute the temperature-dependent factors for the rates.
    real_t t93 = powf(t9, THIRD);
    real_t t1 = 1 / t9;
    real_t t2 = 1 / t93;
    real_t t3 = t93;
    real_t t4 = t9;
    real_t t5 = t93 * t93 * t93 * t93 * t93;
    real_t t6 = logf(t9);

    for (int i = 0; i < number_reactions; i++) {
        Rate[i] =
            Prefac[i] * expf(P0[i] + t1 * P1[i] + t2 * P2[i] + t3 * P3[i] +
                             t4 * P4[i] + t5 * P5[i] + t6 * P6[i]);
    }

    /*
     * Begin the time integration from t=0 to tmax. Rather than t=0 we start at
     * some very small value of t.
     */
    real_t t = 1.0e-16;      // The current integration time
    real_t dt = dt_init;     // The current integration timestep
    real_t prevdt = dt_init; // The integration timestep from the previous step

    real_t dE = 0;

    // Main time integration loop
    while (t < t_max) {
        // Compute the fluxes from the previously-computed rates and the current
        // abundances
        for (int i = 0; i < number_reactions; i++) {
            int nr = NumReactingSpecies[i];
            switch (nr) {
            case 1:
                Flux[i] = Rate[i] * Y[Reactant1[i]];
                break;
            case 2:
                Flux[i] = Rate[i] * Y[Reactant1[i]] * Y[Reactant2[i]];
                break;
            case 3:
                Flux[i] = Rate[i] * Y[Reactant1[i]] * Y[Reactant2[i]] *
                          Y[Reactant3[i]];
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

        // Calculate the change in temperature caused by nuclear burning.
        real_t tmp = 0;
        for (int i = 0; i < number_reactions; i++) {
            tmp += Flux[i] * Q[i];
        }
        tmp *= dt;
        dE += tmp * ECON;

        // Increment the integration time and set the new timestep
        t += dt;
        integration_steps++;
        prevdt = dt;
        dt = compute_timestep(prevdt, t, t_max);
    }

    // Placeholder for proper dT/mC calc.
    // ~ t9 += dE / (rho*volume * c)
    // ...and account for the div by 1e9
    t9 += (dE / (1e05 * 1e9));

    real_t_val[0] = t9;

    return EXIT_SUCCESS;
}

int check_asy(real_t Fminus, real_t Y, real_t dt) {
    if (Y > 0.0f && Fminus * dt / Y > 1.0f) {
        return 1;
    } else {
        return 0;
    }
}

real_t asymptotic_update(real_t Fplus, real_t Fminus, real_t Y, real_t dt) {
    return (Y + Fplus * dt) / (1.0f + Fminus * dt / Y); // Sophia He formula
}

real_t euler_update(real_t FplusSum, real_t FminusSum, real_t dt) {
    return (FplusSum - FminusSum) * dt;
}

real_t compute_timestep(real_t prevdt, real_t t, real_t tmax) {
    real_t dt;
    if (t == 0.0f) {
        dt = 1.0e-20;
    } else {
        dt = 0.1f * t;
    }

    // Prevent final integration step from overstepping tmax
    // if(t+dt > tmax) dt = tmax - t;

    return dt;
}

real_t compute_keff(real_t Fminus, real_t Y) {
    if (Y > 0) {
        return Fminus / Y;
    } else {
        return 0.0f;
    }
}

// Better data organization could massively improve the efficiency of this.
int tnn_data_preprocess(struct tnn**** tnn, struct rt_hydro_mesh* mesh,
                        struct simulation_properties* sim_prop) {
    for (int i = 0; i < sim_prop->resolution[0]; i++) {
        for (int j = 0; j < sim_prop->resolution[1]; j++) {
            for (int k = 0; k < sim_prop->resolution[2]; k++) {
                tnn[i][j][k]->f->rho = mesh->density[i][j][k];
                tnn[i][j][k]->f->t9 = mesh->temp[i][j][k] / 1e9;
                tnn[i][j][k]->f->t_max = mesh->dt;
            }
        }
    }
    return EXIT_SUCCESS;
}

int tnn_data_postprocess(struct tnn**** tnn, struct rt_hydro_mesh* mesh,
                         struct simulation_properties* sim_prop) {
    for (int i = 0; i < sim_prop->resolution[0]; i++) {
        for (int j = 0; j < sim_prop->resolution[1]; j++) {
            for (int k = 0; k < sim_prop->resolution[2]; k++) {
                mesh->temp[i][j][k] = tnn[i][j][k]->f->t9 * 1e9;
            }
        }
    }
    return EXIT_SUCCESS;
}

int problem_parameters_update(struct problem_parameters* params,
                              struct rate_library* rates, struct tnn* network) {

    params->f_plus_total = 0;
    params->f_minus_total = 0;

    int* temp_int1 =
        calloc(network->info->number_species * rates->number_reactions / 2,
               sizeof(int));
    int* temp_int2 =
        calloc(network->info->number_species * rates->number_reactions / 2,
               sizeof(int));

    reaction_mask_update(params->reaction_mask, rates, network, params,
                         temp_int1, temp_int2);

    params->f_plus_isotope_cut[0] = params->f_plus_number[0];
    params->f_minus_isotope_cut[0] = params->f_minus_number[0];
    for (int i = 1; i < network->info->number_species; i++) {
        params->f_plus_isotope_cut[i] =
            params->f_plus_number[i] + params->f_plus_isotope_cut[i - 1];
        params->f_minus_isotope_cut[i] =
            params->f_minus_number[i] + params->f_minus_isotope_cut[i - 1];
    }

    int current_iso = 0;
    for (int i = 0; i < params->f_plus_total; i++) {
        params->f_plus_isotope_idx[i] = current_iso;
        if (i == (params->f_plus_isotope_cut[current_iso] - 1)) {
            current_iso++;
        }
    }

    current_iso = 0;
    for (int i = 0; i < params->f_minus_total; i++) {
        params->f_minus_isotope_idx[i] = current_iso;
        if (i == (params->f_minus_isotope_cut[current_iso] - 1))
            current_iso++;
    }

    for (int i = 0; i < params->f_plus_total; i++) {
        params->f_plus_map[i] = temp_int1[i];
    }

    for (int i = 0; i < params->f_minus_total; i++) {
        params->f_minus_map[i] = temp_int2[i];
    }

    // Populate the params->f_plus_min and params->f_plus_max arrays
    params->f_plus_min[0] = 0;
    params->f_plus_max[0] = params->f_plus_number[0] - 1;
    for (int i = 1; i < network->info->number_species; i++) {
        params->f_plus_min[i] = params->f_plus_max[i - 1] + 1;
        params->f_plus_max[i] =
            params->f_plus_min[i] + params->f_plus_number[i] - 1;
    }
    // Populate the params->f_minus_min and params->f_minus_max arrays
    params->f_minus_min[0] = 0;
    params->f_minus_max[0] = params->f_minus_number[0] - 1;
    for (int i = 1; i < network->info->number_species; i++) {
        params->f_minus_min[i] = params->f_minus_max[i - 1] + 1;
        params->f_minus_max[i] =
            params->f_minus_min[i] + params->f_minus_number[i] - 1;
    }

    int temp_count_plus = 0;
    int temp_count_minus = 0;
    for (int i = 0; i < network->info->number_species; i++) {
        for (int j = 0; j < rates->number_reactions; j++) {
            if (params->reaction_mask[i][j] > 0) {
                params->f_plus_factor[temp_count_plus] =
                    (real_t)params->reaction_mask[i][j];
                temp_count_plus++;
            } else if (params->reaction_mask[i][j] < 0) {
                params->f_minus_factor[temp_count_minus] =
                    -(real_t)params->reaction_mask[i][j];
                temp_count_minus++;
            }
        }
    }

    free(temp_int1);
    free(temp_int2);

    return EXIT_SUCCESS;
}

void reaction_mask_update(int** mask, struct rate_library* rates,
                          struct tnn* network,
                          struct problem_parameters* params, int* temp_int1,
                          int* temp_int2) {

    int increment_plus = 0;
    int increment_minus = 0;

    for (int i = 0; i < network->info->number_species; i++) {
        int total = 0;
        int f_plus_number = 0;
        int f_minus_number = 0;

        for (int j = 0; j < rates->number_reactions; j++) {
            int l_total = 0;
            int r_total = 0;

            for (int k = 0; k < rates->num_react_species[j]; k++) {
                if (network->iptr->z[i] == rates->reactant_z[j][k] &&
                    network->iptr->n[i] == rates->reactant_n[j][k])
                    l_total++;
            }

            for (int k = 0; k < rates->num_products[j]; k++) {
                if (network->iptr->z[i] == rates->product_z[j][k] &&
                    network->iptr->n[i] == rates->product_n[j][k])
                    r_total++;
            }

            total = l_total - r_total;

            if (total > 0) {
                f_minus_number++;
                mask[i][j] = -total;
                temp_int2[increment_minus + f_minus_number - 1] = j;
            } else if (total < 0) { // Contributes to F+ for this isotope
                f_plus_number++;
                mask[i][j] = -total;
                temp_int1[increment_plus + f_plus_number - 1] = j;
            } else {
                mask[i][j] = 0;
            }
        }

        params->f_plus_total += f_plus_number;
        params->f_minus_total += f_minus_number;

        params->f_plus_number[i] = f_plus_number;
        params->f_minus_number[i] = f_minus_number;

        increment_plus += f_plus_number;
        increment_minus += f_minus_number;
    }
}

int tnn_integrate_network(struct problem_parameters* params,
                          struct rate_library* rates, struct tnn* network,
                          struct simulation_properties* sim_prop) {
#ifdef __MP_ROCM
    // Commented out because we want it to run regardless...
    // if (options.rocm_accel) {
    //     printf("THERMONUCLEAR ROCM KERNEL NOT IMPLEMENTED\n");
    //     return EXIT_FAILURE;
    // } else {
    if (tnn_integration_kernel(
            rates->p0, rates->p1, rates->p2, rates->p3, rates->p4, rates->p5,
            rates->p6, rates->prefactor, rates->q_value, rates->rate,
            rates->flux, params->f_plus, params->f_minus, params->f_plus_factor,
            params->f_minus_factor, params->f_plus_sum, params->f_minus_sum,
            params->f_plus_max, params->f_minus_max, params->f_plus_map,
            params->f_minus_map, network->fptr->y, rates->num_react_species,
            rates->reactant_1, rates->reactant_2, rates->reactant_3,
            network->info->number_species, rates->number_reactions,
            params->f_plus_total, params->f_minus_total, network->f->t9,
            network->f->t_max, network->f->dt_init) == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }
    // }
#else
    // tmp? More info is in the neutrino kernel (in the same place).
    real_t* tmp = malloc(1 * sizeof(real_t));
    tmp[0] = network->f->t9;
    if (tnn_integration_kernel(
            tmp, NULL, rates->p0, rates->p1, rates->p2, rates->p3, rates->p4,
            rates->p5, rates->p6, params->prefactor, rates->q_value,
            rates->rate, rates->flux, params->f_plus, params->f_minus,
            params->f_plus_factor, params->f_minus_factor, params->f_plus_sum,
            params->f_minus_sum, params->f_plus_max, params->f_minus_max,
            params->f_plus_map, params->f_minus_map, network->fptr->y,
            rates->num_react_species, rates->reactant_1, rates->reactant_2,
            rates->reactant_3, network->info->number_species,
            rates->number_reactions, params->f_plus_total,
            params->f_minus_total, network->f->t9, network->f->t_max,
            network->f->dt_init) == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }
    network->f->t9 = tmp[0];
    free(tmp);
#endif

    return EXIT_SUCCESS;
}

int tnn_kernel_trigger(struct rate_library* rates, struct tnn**** network,
                       struct problem_parameters* params,
                       struct simulation_properties* sim_prop) {

    for (int i = 0; i < sim_prop->resolution[0]; i++) {
        for (int j = 0; j < sim_prop->resolution[1]; j++) {
            for (int k = 0; k < sim_prop->resolution[2]; k++) {
                // This per-network preprocessing should occur in
                // the compute kernel.
                real_t density[3];
                density[0] = 1.0f;
                density[1] = network[i][j][k]->f->rho;
                density[2] =
                    network[i][j][k]->f->rho * network[i][j][k]->f->rho;
                problem_parameters_update(params, rates, network[i][j][k]);
                for (int n = 0; n < rates->number_reactions; n++) {
                    params->prefactor[n] =
                        rates->prefactor[n] *
                        (density[rates->num_react_species[n] - 1]);
                }
                if (tnn_integrate_network(params, rates, network[i][j][k],
                                          sim_prop) == EXIT_FAILURE) {
                    return EXIT_FAILURE;
                }
            }
        }
    }
    return EXIT_SUCCESS;
}
