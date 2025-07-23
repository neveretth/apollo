#include "kernel.h"

#include <math.h>
#include <stdlib.h>

int hydro_integration_kernel(real_t*** temp, real_t*** density, real_t volume,
                             real_t h, real_t dt, real_t t_end, int* dim) {
    // This really should be passed, not created each time.
    real_t*** delta_temp = malloc(dim[0] * sizeof(real_t*));
    for (int i = 0; i < dim[0]; i++) {
        delta_temp[i] = malloc(dim[1] * sizeof(real_t*));
        for (int j = 0; j < dim[1]; j++) {
            delta_temp[i][j] = malloc(dim[2] * sizeof(real_t));
        }
    }

    // For now just assume it's all even. Volume is for each cell.
    real_t area = volume / dim[2]; // Placeholder area of interaction.

    real_t c = 1e8; // Placeholder contribution value.
    real_t ntd = 0; // Nuclear burning temp diff (assume negligible for now)

    real_t t = 0;

    // Does one integration (single step if t_end = 0)
    do {
        for (int i = 0; i < dim[0]; i++) {
            for (int j = 0; j < dim[1]; j++) {
                for (int k = 0; k < dim[2]; k++) {
                    real_t pdvh = 1 / (density[i][j][k] * volume * c);
                    real_t temp_diff = -(6 * temp[i][j][k]);

                    // As with the linear integration, towards idx=0 is positive
                    temp_diff += temp[i - (i > 0)][j][k];
                    temp_diff += temp[i + (i < dim[0] - 1)][j][k];
                    temp_diff += temp[i][j - (j > 0)][k];
                    temp_diff += temp[i][j + (j < dim[1] - 1)][k];
                    temp_diff += temp[i][j][k - (k > 0)];
                    temp_diff += temp[i][j][k + (k < dim[2] - 1)];

                    delta_temp[i][j][k] =
                        dt * (pdvh * (ntd + (h * area * temp_diff)));
                }
            }
        }
        for (int i = 0; i < dim[0]; i++) {
            for (int j = 0; j < dim[1]; j++) {
                for (int k = 0; k < dim[2]; k++) {
                    temp[i][j][k] += delta_temp[i][j][k];
                }
            }
        }
        t += dt;
    } while (t < t_end);

    for (int i = 0; i < dim[0]; i++) {
        for (int j = 0; j < dim[1]; j++) {
            free(delta_temp[i][j]);
        }
        free(delta_temp[i]);
    }
    free(delta_temp);

    return EXIT_SUCCESS;
}

int hydro_data_preprocess() { return EXIT_SUCCESS; }

int hydro_integrate_mesh(struct simulation_properties sim_prop,
                         struct rt_hydro_mesh* mesh,
                         struct option_values options) {
#ifdef __MP_ROCM
    if (hydro_integration_kernel(mesh->temp, mesh->density, mesh->volume,
                                 mesh->h, mesh->dt, mesh->t_end,
                                 mesh->dim) == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }
    // There is no need for a hydro ROCM kernel (right now)
    // if (options.rocm_accel) {
    //     printf("HYDRO ROCM KERNEL NOT IMPLEMENTED\n");
    // } else {
    // }
#else
    if (hydro_integration_kernel(mesh->temp, mesh->density, mesh->volume,
                                 mesh->h, mesh->dt, mesh->t_end,
                                 mesh->dim) == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }
#endif
    return EXIT_SUCCESS;
}
