#include "hydro-types.h"
#include "../types.h"

#include <stdio.h>
#include <stdlib.h>

int hydro_mesh_destroy(struct hydro_mesh** __src) {
    struct hydro_mesh* src = *__src;
    free(src->temp);
    free(src->density);
    free(src);
    *__src = NULL;
    return EXIT_SUCCESS;
}

void hydro_mesh_print(struct hydro_mesh* mesh) {
    printf("Temp: ");
    for (int i = 0; i < mesh->dim; i++) {
        printf("%f ", mesh->temp[i]);
    }
    printf("\n");
    printf("Density: ");
    for (int i = 0; i < mesh->dim; i++) {
        printf("%f ", mesh->density[i]);
    }
    printf("\n");
}

void flat_hydro_mesh_print(struct flat_hydro_mesh* mesh) {
    printf("Temp grid: \n");
    for (int i = 0; i < mesh->dim[0]; i++) {
        for (int j = 0; j < mesh->dim[1]; j++) {
            printf(" %8.3f ", mesh->temp[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int flat_hydro_mesh_destroy(struct flat_hydro_mesh** __src) {
    struct flat_hydro_mesh* src = *__src;
    freenptr((void*)src->temp, src->dim[0]);
    freenptr((void*)src->density, src->dim[0]);
    free(src);
    *__src = NULL;
    return EXIT_SUCCESS;
}

int rt_hydro_mesh_destroy(struct rt_hydro_mesh** __src) {
    struct rt_hydro_mesh* src = *__src;
    for (int i = 0; i < src->dim[0]; i++) {
        for (int j = 0; j < src->dim[1]; j++) {
            free(src->temp[i][j]);
            free(src->density[i][j]);
            free(src->entropy[i][j]);
        }
        free(src->temp[i]);
        free(src->density[i]);
        free(src->entropy[i]);
    }
    free(src->temp);
    free(src->density);
    free(src);
    *__src = NULL;
    return EXIT_SUCCESS;
}
