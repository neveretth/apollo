#include "kernel-driver.h"

#include "hip-util.h"
#include "kernel.h"

#include <stdlib.h>

int integrate_network(struct rate_library* rates, struct thermo_network* network,
                      struct problem_parameters* params,
                      struct option_values options) {

    if (options.rocm_accel) {
        struct hipDeviceProp_t* device = get_hip_device();

        int error;
        struct dim3 blockdim = {256, 1, 1};
        struct dim3 griddim = {1, 1, 1};
        int sharedmem_allocation = 0 * sizeof(float);

        void* d_p0;
        hipMalloc(&d_p0, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_p0, rates->p0, sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_p1;
        hipMalloc(&d_p1, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_p1, rates->p1, sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_p2;
        hipMalloc(&d_p2, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_p2, rates->p2, sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_p3;
        hipMalloc(&d_p3, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_p3, rates->p3, sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_p4;
        hipMalloc(&d_p4, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_p4, rates->p4, sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_p5;
        hipMalloc(&d_p5, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_p5, rates->p5, sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_p6;
        hipMalloc(&d_p6, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_p6, rates->p6, sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_prefactor;
        hipMalloc(&d_prefactor, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_prefactor, rates->prefactor,
                  sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_q_value;
        hipMalloc(&d_q_value, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_q_value, rates->q_value,
                  sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_rate;
        hipMalloc(&d_rate, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_rate, rates->rate, sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_flux;
        hipMalloc(&d_flux, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_flux, rates->flux, sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_f_plus;
        hipMalloc(&d_f_plus, sizeof(float) * params->f_plus_total);
        hipMemcpy(d_f_plus, params->f_plus,
                  sizeof(float) * params->f_plus_total, hipMemcpyHostToDevice);
        void* d_f_minus;
        hipMalloc(&d_f_minus, sizeof(float) * params->f_minus_total);
        hipMemcpy(d_f_minus, params->f_minus,
                  sizeof(float) * params->f_plus_total, hipMemcpyHostToDevice);
        void* d_f_plus_factor;
        hipMalloc(&d_f_plus_factor, sizeof(float) * params->f_minus_total);
        hipMemcpy(d_f_plus_factor, params->f_plus_factor,
                  sizeof(float) * params->f_plus_total, hipMemcpyHostToDevice);
        void* d_f_minus_factor;
        hipMalloc(&d_f_minus_factor, sizeof(float) * params->f_minus_total);
        hipMemcpy(d_f_minus_factor, params->f_minus_factor,
                  sizeof(float) * params->f_minus_total, hipMemcpyHostToDevice);
        void* d_f_plus_sum;
        hipMalloc(&d_f_plus_sum, sizeof(float) * network->info->number_species);
        hipMemcpy(d_f_plus_sum, params->f_plus_sum,
                  sizeof(float) * network->info->number_species,
                  hipMemcpyHostToDevice);
        void* d_f_minus_sum;
        hipMalloc(&d_f_minus_sum, sizeof(float) * network->info->number_species);
        hipMemcpy(d_f_minus_sum, params->f_minus_sum,
                  sizeof(float) * network->info->number_species,
                  hipMemcpyHostToDevice);
        void* d_f_plus_max;
        hipMalloc(&d_f_plus_max, sizeof(float) * network->info->number_species);
        hipMemcpy(d_f_plus_max, params->f_plus_max,
                  sizeof(float) * network->info->number_species,
                  hipMemcpyHostToDevice);
        void* d_f_minus_max;
        hipMalloc(&d_f_minus_max, sizeof(float) * network->info->number_species);
        hipMemcpy(d_f_minus_max, params->f_minus_max,
                  sizeof(float) * network->info->number_species,
                  hipMemcpyHostToDevice);
        void* d_f_plus_map;
        hipMalloc(&d_f_plus_map, sizeof(float) * params->f_plus_total);
        hipMemcpy(d_f_plus_map, params->f_plus_map,
                  sizeof(float) * params->f_plus_total, hipMemcpyHostToDevice);
        void* d_f_minus_map;
        hipMalloc(&d_f_minus_map, sizeof(float) * params->f_minus_total);
        hipMemcpy(d_f_minus_map, params->f_minus_map,
                  sizeof(float) * params->f_minus_total, hipMemcpyHostToDevice);
        void* d_y;
        hipMalloc(&d_y, sizeof(float) * network->info->number_species);
        hipMemcpy(d_y, network->fptr->y, sizeof(float) * network->info->number_species,
                  hipMemcpyHostToDevice);
        void* d_num_react_species;
        hipMalloc(&d_num_react_species,
                  sizeof(float) * rates->number_reactions);
        hipMemcpy(d_num_react_species, rates->num_react_species,
                  sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_reactant_1;
        hipMalloc(&d_reactant_1, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_reactant_1, rates->reactant_1,
                  sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_reactant_2;
        hipMalloc(&d_reactant_2, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_reactant_2, rates->reactant_2,
                  sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_reactant_3;
        hipMalloc(&d_reactant_3, sizeof(float) * rates->number_reactions);
        hipMemcpy(d_reactant_3, rates->reactant_3,
                  sizeof(float) * rates->number_reactions,
                  hipMemcpyHostToDevice);
        void* d_int_val;
        int* int_val = malloc(sizeof(int) * 5);
        int_val[0] = network->info->number_species;
        int_val[1] = rates->number_reactions;
        int_val[2] = params->f_plus_total;
        int_val[3] = params->f_minus_total;
        int_val[4] = options.halt;
        hipMalloc(&d_int_val, sizeof(int) * 5);
        hipMemcpy(d_int_val, int_val, sizeof(float) * 5, hipMemcpyHostToDevice);
        void* d_float_val;
        float* float_val = malloc(sizeof(float) * 3);
        float_val[0] = network->f->t9;
        float_val[1] = network->f->t_max;
        float_val[2] = network->f->dt_init;
        hipMalloc(&d_float_val, sizeof(float) * 3);
        hipMemcpy(d_float_val, float_val, sizeof(float) * 3,
                  hipMemcpyHostToDevice);

        free(int_val);
        free(float_val);

        // This is particularly unpleasant.
        // There's _supposedly_ a way to do this with hipLaunchKernelGGL(),
        // but alas HIP documentation is like fondant:
        // pretty, but it ruins your day.
        void** args = malloc(sizeof(void*) * 28);
        args[0] = &d_p0;
        args[1] = &d_p1;
        args[2] = &d_p2;
        args[3] = &d_p3;
        args[4] = &d_p4;
        args[5] = &d_p5;
        args[6] = &d_p6;
        args[7] = &d_prefactor;
        args[8] = &d_q_value;
        args[9] = &d_rate;
        args[10] = &d_flux;
        args[11] = &d_f_plus;
        args[12] = &d_f_minus;
        args[13] = &d_f_plus_factor;
        args[14] = &d_f_minus_factor;
        args[15] = &d_f_plus_sum;
        args[16] = &d_f_minus_sum;
        args[17] = &d_f_plus_max;
        args[18] = &d_f_minus_max;
        args[19] = &d_f_plus_map;
        args[20] = &d_f_minus_map;
        args[21] = &d_y;
        args[22] = &d_num_react_species;
        args[23] = &d_reactant_1;
        args[24] = &d_reactant_2;
        args[25] = &d_reactant_3;
        args[26] = &d_int_val;
        args[27] = &d_float_val;

        if ((error = hipLaunchKernel(integration_kernel_hip, griddim, blockdim,
                                     args, sharedmem_allocation,
                                     hipStreamDefault)) != hipSuccess) {
            printf("==apollo== encountered error launching kernel: %i\n",
                   error);
            return EXIT_FAILURE;
        }
        hipMemcpy(rates->rate, d_rate, sizeof(float) * rates->number_reactions,
                  hipMemcpyDeviceToHost);
        hipMemcpy(rates->flux, d_flux, sizeof(float) * rates->number_reactions,
                  hipMemcpyDeviceToHost);
        hipMemcpy(params->f_plus, d_f_plus,
                  sizeof(float) * params->f_plus_total, hipMemcpyDeviceToHost);
        hipMemcpy(params->f_minus, d_f_minus,
                  sizeof(float) * params->f_plus_total, hipMemcpyDeviceToHost);
        hipMemcpy(params->f_plus_factor, d_f_plus_factor,
                  sizeof(float) * params->f_plus_total, hipMemcpyDeviceToHost);
        hipMemcpy(params->f_minus_factor, d_f_minus_factor,
                  sizeof(float) * params->f_minus_total, hipMemcpyDeviceToHost);
        hipMemcpy(params->f_plus_sum, d_f_plus_sum,
                  sizeof(float) * network->info->number_species,
                  hipMemcpyDeviceToHost);
        hipMemcpy(params->f_minus_sum, d_f_minus_sum,
                  sizeof(float) * network->info->number_species,
                  hipMemcpyDeviceToHost);
        hipMemcpy(params->f_plus_max, d_f_plus_max,
                  sizeof(float) * network->info->number_species,
                  hipMemcpyDeviceToHost);
        hipMemcpy(params->f_minus_max, d_f_minus_max,
                  sizeof(float) * network->info->number_species,
                  hipMemcpyDeviceToHost);
        hipMemcpy(params->f_plus_map, d_f_plus_map,
                  sizeof(float) * params->f_plus_total, hipMemcpyDeviceToHost);
        hipMemcpy(params->f_minus_map, d_f_minus_map,
                  sizeof(float) * params->f_minus_total, hipMemcpyDeviceToHost);
        hipMemcpy(network->fptr->y, d_y, sizeof(float) * network->info->number_species,
                  hipMemcpyDeviceToHost);
        for (int i = 0; i < 34; i++) {
            // I'm just that goddamn lazy.
            hipFree(args[i]);
        }
        free(args);
        free(device);
    } else {
        if (integration_kernel(
                rates->p0, rates->p1, rates->p2, rates->p3, rates->p4,
                rates->p5, rates->p6, rates->prefactor, rates->q_value,
                rates->rate, rates->flux, params->f_plus, params->f_minus,
                params->f_plus_factor, params->f_minus_factor,
                params->f_plus_sum, params->f_minus_sum, params->f_plus_max,
                params->f_minus_max, params->f_plus_map, params->f_minus_map,
                network->fptr->y, rates->num_react_species, rates->reactant_1,
                rates->reactant_2, rates->reactant_3, network->info->number_species,
                rates->number_reactions, params->f_plus_total,
                params->f_minus_total, network->f->t9, network->f->t_max,
                network->f->dt_init, options.halt) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
