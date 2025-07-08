#include "kernel-driver.h"

#include "../hip-util.h"
#include "../kernel/kernel.h"

#include <stdlib.h>

int tnn_integrate_network(struct rate_library* rates, struct tnn* network,
                          struct problem_parameters* params,
                          struct option_values options) {

    if (options.rocm_accel) {
        printf("==apollo== You are running the NON-rocm apollo, try apollo-rocm!\n");
        return EXIT_FAILURE;
    } else {
        if (tnn_integration_kernel(
                rates->p0, rates->p1, rates->p2, rates->p3, rates->p4,
                rates->p5, rates->p6, rates->prefactor, rates->q_value,
                rates->rate, rates->flux, params->f_plus, params->f_minus,
                params->f_plus_factor, params->f_minus_factor,
                params->f_plus_sum, params->f_minus_sum, params->f_plus_max,
                params->f_minus_max, params->f_plus_map, params->f_minus_map,
                network->fptr->y, rates->num_react_species, rates->reactant_1,
                rates->reactant_2, rates->reactant_3,
                network->info->number_species, rates->number_reactions,
                params->f_plus_total, params->f_minus_total, network->f->t9,
                network->f->t_max, network->f->dt_init,
                options.halt) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int neunet_integrate_network(struct neunet* network,
                             struct option_values options) {
    if (options.rocm_accel) {
        printf("==apollo== You are running the NON-rocm apollo, try apollo-rocm!\n");
        return EXIT_FAILURE;
    } else {
        if (neunet_integration_kernel(
                network->info->rate_in, network->info->rate_out,
                network->fptr->n_old, network->fptr->ec, network->fptr->dv,
                network->f->dt, network->f->t_end, network->f->EpsA,
                network->f->EpsR, network->f->g_a, network->f->g_b,
                network->f->g_c, network->info->num_groups,
                options.halt) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int hydro_integrate_mesh(struct hydro_mesh* mesh,
                         struct option_values options) {
    if (options.rocm_accel) {
        printf("==apollo== You are running the NON-rocm apollo, try apollo-rocm!\n");
        return EXIT_FAILURE;
    } else {
        if (hydro_integration_kernel(mesh->temp, mesh->density, mesh->volume,
                                     mesh->h, mesh->dt, mesh->t_end,
                                     mesh->dim) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
