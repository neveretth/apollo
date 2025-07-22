#include "../kernel-driver.h"

#include "../../hip-util.h"
#include "../../kernel/thermonuclear/kernel.h"
#include "../../kernel/neutrino/kernel.h"
#include "../../kernel/hydro/kernel.h"

#include <stdlib.h>

int tnn_integrate_network(struct rate_library* rates, struct tnn* network,
                          struct problem_parameters* params,
                          struct option_values options) {

    if (options.rocm_accel) {
        struct hipDeviceProp_t* device = get_hip_device();

        int error;
        struct dim3 blockdim = {256, 1, 1};
        struct dim3 griddim = {1, 1, 1};
        int sharedmem_allocation = 0 * sizeof(float);

        int number_args = 28;
        void** args = malloc(sizeof(void*) * number_args);
        for (int i = 0; i < number_args; i++) {
            args[i] = malloc(sizeof(void*));
        }
        devbuf_create(args[0], rates->p0,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[1], rates->p1,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[2], rates->p2,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[3], rates->p3,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[4], rates->p4,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[5], rates->p5,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[6], rates->p6,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[7], rates->prefactor,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[8], rates->q_value,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[9], rates->rate,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[10], rates->flux,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[11], params->f_plus,
                      params->f_plus_total * sizeof(float));
        devbuf_create(args[12], params->f_minus,
                      params->f_minus_total * sizeof(float));
        devbuf_create(args[13], params->f_plus_factor,
                      params->f_plus_total * sizeof(float));
        devbuf_create(args[14], params->f_minus_factor,
                      params->f_minus_total * sizeof(float));
        devbuf_create(args[15], params->f_plus_sum,
                      network->info->number_species * sizeof(float));
        devbuf_create(args[16], params->f_minus_sum,
                      network->info->number_species * sizeof(float));
        devbuf_create(args[17], params->f_plus_max,
                      network->info->number_species * sizeof(float));
        devbuf_create(args[18], params->f_minus_max,
                      network->info->number_species * sizeof(float));
        devbuf_create(args[19], params->f_plus_map,
                      params->f_plus_total * sizeof(float));
        devbuf_create(args[20], params->f_minus_map,
                      params->f_minus_total * sizeof(float));
        devbuf_create(args[21], network->fptr->y,
                      network->info->number_species * sizeof(float));
        devbuf_create(args[22], rates->num_react_species,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[23], rates->reactant_1,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[24], rates->reactant_2,
                      rates->number_reactions * sizeof(float));
        devbuf_create(args[25], rates->reactant_3,
                      rates->number_reactions * sizeof(float));

        void* d_int_val;
        int* int_val = malloc(sizeof(int) * 5);
        int_val[0] = network->info->number_species;
        int_val[1] = rates->number_reactions;
        int_val[2] = params->f_plus_total;
        int_val[3] = params->f_minus_total;
        int_val[4] = options.halt;
        devbuf_create(args[26], int_val, 5 * sizeof(int));
        free(int_val);
        
        float* float_val = malloc(sizeof(float) * 3);
        float_val[0] = network->f->t9;
        float_val[1] = network->f->t_max;
        float_val[2] = network->f->dt_init;
        devbuf_create(args[27], float_val, 3 * sizeof(float));
        free(float_val);

        if ((error = hipLaunchKernel(integration_kernel_hip, griddim, blockdim,
                                     args, sharedmem_allocation,
                                     hipStreamDefault)) != hipSuccess) {
            printf("==apollo== encountered error launching kernel: %i\n",
                   error);
            return EXIT_FAILURE;
        }
        devbuf_read(rates->rate, args[9],
                    sizeof(float) * rates->number_reactions);
        devbuf_read(rates->flux, args[10],
                    sizeof(float) * rates->number_reactions);
        devbuf_read(params->f_plus, args[11],
                    sizeof(float) * params->f_plus_total);
        devbuf_read(params->f_minus, args[12],
                    sizeof(float) * params->f_plus_total);
        devbuf_read(params->f_plus_factor, args[13],
                    sizeof(float) * params->f_plus_total);
        devbuf_read(params->f_minus_factor, args[14],
                    sizeof(float) * params->f_minus_total);
        devbuf_read(params->f_plus_sum, args[15],
                    sizeof(float) * network->info->number_species);
        devbuf_read(params->f_minus_sum, args[16],
                    sizeof(float) * network->info->number_species);
        devbuf_read(params->f_plus_max, args[17],
                    sizeof(float) * network->info->number_species);
        devbuf_read(params->f_minus_max, args[18],
                    sizeof(float) * network->info->number_species);
        devbuf_read(params->f_plus_map, args[19],
                    sizeof(float) * params->f_plus_total);
        devbuf_read(params->f_minus_map, args[20],
                    sizeof(float) * params->f_minus_total);
        devbuf_read(network->fptr->y, args[21],
                    sizeof(float) * network->info->number_species);
        for (int i = 0; i < 34; i++) {
            hipFree(args[i]);
        }
        free(args);
        free(device);
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
        printf("Not implemented\n");
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
        printf("Not implemented\n");
    } else {
        if (hydro_integration_kernel(mesh->temp, mesh->density, mesh->volume,
                                     mesh->h, mesh->dt, mesh->t_end,
                                     mesh->dim) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int hydro_integrate_flat_mesh(struct flat_hydro_mesh* mesh,
                         struct option_values options) {
    if (options.rocm_accel) {
        printf("Not implemented\n");
    } else {
        if (flat_hydro_integration_kernel(mesh->temp, mesh->density,
                                     mesh->volume, mesh->h, mesh->dt,
                                     mesh->t_end, mesh->dim) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

int hydro_integrate_rt_mesh(struct rt_hydro_mesh* mesh,
                            struct option_values options) {
    if (options.rocm_accel) {
        printf("==apollo== You are running the NON-rocm apollo, try "
               "apollo-rocm!\n");
    } else {
        if (rt_hydro_integration_kernel(mesh->temp, mesh->density, mesh->volume,
                                        mesh->h, mesh->dt, mesh->t_end,
                                        mesh->dim) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
