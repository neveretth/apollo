#include "kernel.h"

#include <math.h>
#include <stdlib.h>

int advout_entropy(struct advout_t* data, struct rt_hydro_mesh* mesh,
                   struct simulation_properties* sim_prop) {
    for (int k = 0; k < mesh->dim[2]; k++) {
        for (int j = 0; j < mesh->dim[1]; j++) {
            for (int i = 0; i < mesh->dim[0]; i++) {
                data->entropy[k][j][i] +=
                    1e8 * 2 *
                    ((mesh->temp[k][j][i] - data->prev_temp[k][j][i]) /
                     (mesh->temp[k][j][i] + data->prev_temp[k][j][i]));
                data->prev_temp[k][j][i] = mesh->temp[k][j][i];
            }
        }
    }
    return EXIT_SUCCESS;
}

struct advout_t* advout_data_create(struct rt_hydro_mesh* mesh,
                                    struct simulation_properties sim_prop) {
    struct advout_t* data = malloc(sizeof(struct advout_t));

    int flat_len = mesh->dim[XDIM] * mesh->dim[YDIM] * mesh->dim[ZDIM];

    data->entropy_block = malloc(flat_len * sizeof(real_t));
    data->prev_temp_block = malloc(flat_len * sizeof(real_t));

    data->entropy = malloc(mesh->dim[0] * sizeof(real_t*));
    data->prev_temp = malloc(mesh->dim[0] * sizeof(real_t*));
    for (int k = 0; k < mesh->dim[ZDIM]; k++) {
        data->entropy[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
        data->prev_temp[k] = malloc(mesh->dim[YDIM] * sizeof(real_t*));
        for (int j = 0; j < mesh->dim[YDIM]; j++) {
            data->entropy[k][j] = (real_t*)data->entropy_block +
                                  (mesh->dim[XDIM] * j) +
                                  (mesh->dim[XDIM] * mesh->dim[YDIM] * k);
            data->prev_temp[k][j] = (real_t*)data->prev_temp_block +
                                    (mesh->dim[XDIM] * j) +
                                    (mesh->dim[XDIM] * mesh->dim[YDIM] * k);
        }
    }
    return data;
}

int advout_data_setup(struct advout_t* data, struct rt_hydro_mesh* mesh,
                      struct simulation_properties sim_prop) {
    for (int k = 0; k < mesh->dim[2]; k++) {
        for (int j = 0; j < mesh->dim[1]; j++) {
            for (int i = 0; i < mesh->dim[0]; i++) {
                data->entropy[k][j][i] = 0;
                data->prev_temp[k][j][i] = mesh->temp[k][j][i];
            }
        }
    }
    return EXIT_SUCCESS;
}

int advout_data_destroy(struct advout_t** data,
                        struct simulation_properties sim_prop) {
    struct advout_t* src = *data;
    for (int k = 0; k < sim_prop.resolution[ZDIM]; k++) {
        free(src->entropy[k]);
        free(src->prev_temp[k]);
    }
    free(src->entropy);
    free(src->prev_temp);

    free(src->entropy_block);
    free(src->prev_temp_block);

    free(src);
    *data = NULL;
    return EXIT_SUCCESS;
}
