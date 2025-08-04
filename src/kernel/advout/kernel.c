#include "kernel.h"

#include <math.h>
#include <stdlib.h>

int advout_entropy(struct advout_t* data, struct rt_hydro_mesh* mesh,
                   struct simulation_properties* sim_prop) {
    for (int i = 0; i < mesh->dim[0]; i++) {
        for (int j = 0; j < mesh->dim[1]; j++) {
            for (int k = 0; k < mesh->dim[2]; k++) {
                data->entropy[i][j][k] +=
                    1e8 * 2 *
                    ((mesh->temp[i][j][k] - data->prev_temp[i][j][k]) /
                     (mesh->temp[i][j][k] + data->prev_temp[i][j][k]));
                data->prev_temp[i][j][k] = mesh->temp[i][j][k];
            }
        }
    }
    return EXIT_SUCCESS;
}

struct advout_t* advout_data_create(struct rt_hydro_mesh* mesh,
                                    struct simulation_properties sim_prop) {
    struct advout_t* data = malloc(sizeof(struct advout_t));
    data->entropy = malloc(mesh->dim[0] * sizeof(real_t*));
    data->prev_temp = malloc(mesh->dim[0] * sizeof(real_t*));
    for (int i = 0; i < mesh->dim[0]; i++) {
        data->entropy[i] = malloc(mesh->dim[1] * sizeof(real_t*));
        data->prev_temp[i] = malloc(mesh->dim[1] * sizeof(real_t*));
        for (int j = 0; j < mesh->dim[1]; j++) {
            data->entropy[i][j] = malloc(mesh->dim[2] * sizeof(real_t));
            data->prev_temp[i][j] = malloc(mesh->dim[2] * sizeof(real_t));
        }
    }
    return data;
}

int advout_data_setup(struct advout_t* data, struct rt_hydro_mesh* mesh,
                      struct simulation_properties sim_prop) {
    for (int i = 0; i < mesh->dim[0]; i++) {
        for (int j = 0; j < mesh->dim[1]; j++) {
            for (int k = 0; k < mesh->dim[2]; k++) {
                data->entropy[i][j][k] = 0;
                data->prev_temp[i][j][k] = mesh->temp[i][j][k];
            }
        }
    }
    return EXIT_SUCCESS;
}

int advout_data_destroy(struct advout_t** data,
                        struct simulation_properties sim_prop) {
    struct advout_t* src = *data;
    for (int i = 0; i < sim_prop.resolution[0]; i++) {
        for (int j = 0; j < sim_prop.resolution[1]; j++) {
            free(src->entropy[i][j]);
            free(src->prev_temp[i][j]);
        }
        free(src->entropy[i]);
        free(src->prev_temp[i]);
    }
    free(src->entropy);
    free(src->prev_temp);
    free(src);
    *data = NULL;
    return EXIT_SUCCESS;
}
