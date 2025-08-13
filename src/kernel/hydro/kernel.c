#include "kernel.h"

#include <math.h>
#include <stdlib.h>

int hydro_integration_kernel(real_t*** temp, real_t*** density,
                             real_t*** delta_temp, real_t*** pressure,
                             real_t*** velocity, real_t*** mean_mol_mass,
                             real_t volume, real_t h, real_t dt, real_t t_end,
                             int* dim) {

    // For now just assume it's all even. Volume is for each cell.
    real_t area = volume / dim[2]; // Placeholder area of interaction.

    real_t c = 1e8; // Placeholder thermal coefficient.
    real_t ntd = 0; // Nuclear burning temp diff (assume negligible for now)

    real_t t = 0;

    // Does one integration (single step if t_end = 0)
    do {
        for (int k = 0; k < dim[ZDIM]; k++) {
            for (int j = 0; j < dim[YDIM]; j++) {
                for (int i = 0; i < dim[XDIM]; i++) {
                    // TEMPERATURE / HEAT TRANSFER
                    real_t pdvh = 1 / (density[k][j][i] * volume * c);
                    real_t temp_diff = -(6 * temp[k][j][i]);

                    // Towards idx=0 is positive
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
        for (int k = 0; k < dim[ZDIM]; k++) {
            for (int j = 0; j < dim[YDIM]; j++) {
                for (int i = 0; i < dim[XDIM]; i++) {
                    // TEMPERATURE / HEAT TRANSFER
                    temp[k][j][i] += delta_temp[k][j][i];
                    
                    // MEAN MOLECULAR MASS
                    mean_mol_mass[k][j][i] = 1;

                    // PRESSURE
                    pressure[k][j][i] = 1;

                    // VELOCITY
                    velocity[k][j][i + XDIM] = 1;
                    if (dim[YDIM] > 1) {
                        velocity[k][j][i + YDIM] = 1;
                    }
                    if (dim[ZDIM] > 1) {
                        velocity[k][j][i + ZDIM] = 1;
                    }

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
    if (hydro_integration_kernel(
            mesh->temp, mesh->density, mesh->delta_temp, mesh->pressure,
            mesh->velocity, mesh->mean_mol_mass, mesh->volume, mesh->h,
            mesh->dt, mesh->t_end, mesh->dim) == EXIT_FAILURE) {
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
    mesh->pressure_block = malloc(flat_len * sizeof(real_t));
    mesh->mean_mol_mass_block = malloc(flat_len * sizeof(real_t));

    int dimensions =
        (mesh->dim[XDIM] > 1) + (mesh->dim[YDIM] > 1) + (mesh->dim[ZDIM] > 1);
    mesh->velocity_block = malloc(flat_len * dimensions * sizeof(real_t));

    mesh->temp = malloc(mesh->dim[ZDIM] * sizeof(real_t*));
    mesh->density = malloc(mesh->dim[ZDIM] * sizeof(real_t*));
    mesh->delta_temp = malloc(mesh->dim[XDIM] * sizeof(real_t*));
    mesh->pressure = malloc(mesh->dim[XDIM] * sizeof(real_t*));
    mesh->velocity = malloc(mesh->dim[XDIM] * sizeof(real_t*));
    mesh->mean_mol_mass = malloc(mesh->dim[XDIM] * sizeof(real_t*));
    for (int k = 0; k < mesh->dim[ZDIM]; k++) {
        mesh->temp[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
        mesh->density[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
        mesh->delta_temp[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
        mesh->pressure[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
        mesh->velocity[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
        mesh->mean_mol_mass[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
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
            mesh->pressure[k][j] = (real_t*)mesh->pressure_block +
                                   (mesh->dim[XDIM] * j) +
                                   (mesh->dim[XDIM] * mesh->dim[YDIM] * k);
            mesh->mean_mol_mass[k][j] = (real_t*)mesh->mean_mol_mass_block +
                                        (mesh->dim[XDIM] * j) +
                                        (mesh->dim[XDIM] * mesh->dim[YDIM] * k);
            mesh->velocity[k][j] =
                (real_t*)mesh->velocity_block +
                (dimensions * ((mesh->dim[XDIM] * j) +
                               (mesh->dim[XDIM] * mesh->dim[YDIM] * k)));
        }
    }

    return mesh;
}
