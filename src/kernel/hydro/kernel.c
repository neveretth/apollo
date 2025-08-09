#include "kernel.h"

#include <math.h>
#include <stdlib.h>

int hydro_integration_kernel(real_t*** temp, real_t*** density,
                             real_t*** delta_temp, real_t volume, real_t h,
                             real_t dt, real_t t_end, int* dim) {

    // For now just assume it's all even. Volume is for each cell.
    real_t area = volume / dim[2]; // Placeholder area of interaction.

    real_t c = 1e8; // Placeholder thermal coefficient.
    real_t ntd = 0; // Nuclear burning temp diff (assume negligible for now)

    real_t t = 0;

    // Does one integration (single step if t_end = 0)
    do {
        for (int k = 0; k < dim[2]; k++) {
            for (int j = 0; j < dim[1]; j++) {
                for (int i = 0; i < dim[0]; i++) {
                    real_t pdvh = 1 / (density[k][j][i] * volume * c);
                    real_t temp_diff = -(6 * temp[k][j][i]);

                    // As with the linear integration, towards idx=0 is positive
                    temp_diff += temp[k][j][i - (i > 0)];
                    temp_diff += temp[k][j][i + (i < dim[0] - 1)];
                    temp_diff += temp[k][j - (j > 0)][i];
                    temp_diff += temp[k][j + (j < dim[1] - 1)][i];
                    temp_diff += temp[k - (k > 0)][j][i];
                    temp_diff += temp[k + (k < dim[2] - 1)][j][i];

                    delta_temp[k][j][i] =
                        dt * (pdvh * (ntd + (h * area * temp_diff)));
                }
            }
        }
        for (int k = 0; k < dim[2]; k++) {
            for (int j = 0; j < dim[1]; j++) {
                for (int i = 0; i < dim[0]; i++) {
                    temp[k][j][i] += delta_temp[k][j][i];
                }
            }
        }
        t += dt;
    } while (t < t_end);

    return EXIT_SUCCESS;
}

int hydro_data_preprocess() { return EXIT_SUCCESS; }

int hydro_data_postprocess() { return EXIT_SUCCESS; }

int hydro_integrate_mesh(struct rt_hydro_mesh* mesh,
                         struct simulation_properties* sim_prop) {
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
    if (hydro_integration_kernel(mesh->temp, mesh->density, mesh->delta_temp,
                                 mesh->volume, mesh->h, mesh->dt, mesh->t_end,
                                 mesh->dim) == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }
#endif
    return EXIT_SUCCESS;
}

struct rt_hydro_mesh* hydro_mesh_create(struct simulation_properties sim_prop) {
    struct rt_hydro_mesh* mesh = malloc(sizeof(struct rt_hydro_mesh));
    mesh->dim[XDIM] = sim_prop.resolution[XDIM];
    mesh->dim[YDIM] = sim_prop.resolution[YDIM];
    mesh->dim[ZDIM] = sim_prop.resolution[ZDIM];

    int flat_len = mesh->dim[XDIM] * mesh->dim[YDIM] * mesh->dim[ZDIM];

    mesh->temp_block = malloc(flat_len * sizeof(real_t));
    mesh->density_block = malloc(flat_len * sizeof(real_t));
    mesh->delta_temp_block = malloc(flat_len * sizeof(real_t));

    mesh->temp = malloc(mesh->dim[XDIM] * sizeof(real_t*));
    mesh->density = malloc(mesh->dim[XDIM] * sizeof(real_t*));
    mesh->delta_temp = malloc(mesh->dim[XDIM] * sizeof(real_t*));

    mesh->temp = malloc(mesh->dim[ZDIM] * sizeof(real_t*));
    mesh->density = malloc(mesh->dim[ZDIM] * sizeof(real_t*));
    mesh->delta_temp = malloc(mesh->dim[XDIM] * sizeof(real_t*));
    for (int k = 0; k < mesh->dim[ZDIM]; k++) {
        mesh->temp[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
        mesh->density[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
        mesh->delta_temp[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
        for (int j = 0; j < mesh->dim[YDIM]; j++) {
            mesh->temp[k][j] = (real_t*)mesh->temp_block +
                               (mesh->dim[XDIM] * j) +
                               (mesh->dim[XDIM] * mesh->dim[YDIM] * k);
            mesh->density[k][j] = (real_t*)mesh->density_block +
                                  (mesh->dim[XDIM] * j) +
                                  (mesh->dim[XDIM] * mesh->dim[YDIM] * k);
            mesh->delta_temp[k][j] = (real_t*)mesh->delta_temp_block +
                                  (mesh->dim[XDIM] * j) +
                                  (mesh->dim[XDIM] * mesh->dim[YDIM] * k);
        }
    }

    return mesh;
}
