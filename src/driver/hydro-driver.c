#include "hydro-driver.h"
#include "kernel-driver.h"

#include <stdlib.h>

// Debug code for linear_hydro_mesh.
static int linear_hydro_debug(struct option_values options) {
    struct hydro_mesh* mesh = malloc(sizeof(struct hydro_mesh));
    // Represents 100 linear zones.
    mesh->dim = 4;
    mesh->dt = 1e-8;
    mesh->h = 10e+12;
    mesh->volume = 1;
    mesh->temp = malloc(mesh->dim * sizeof(float));
    mesh->density = malloc(mesh->dim * sizeof(float));

    mesh->t_end = 1e-1;
    // mesh->t_end = 1;

    for (int i = 0; i < mesh->dim; i++) {
        mesh->temp[i] = 8e+09 + (rand() % (long)1e08);
        mesh->density[i] = 1; // Even density for debugging.
    }

    hydro_mesh_print(mesh);

    if (hydro_integrate_mesh(mesh, options) == EXIT_FAILURE) {
        hydro_mesh_destroy(&mesh);
        return EXIT_FAILURE;
    }
    hydro_mesh_print(mesh);

    hydro_mesh_destroy(&mesh);
    return EXIT_SUCCESS;
}

// Debug code for flat_hydro_mesh.
static int flat_hydro_debug(struct option_values options) {
    struct flat_hydro_mesh* mesh = malloc(sizeof(struct flat_hydro_mesh));
    mesh->dim[0] = 16;
    mesh->dim[1] = 8;

    mesh->dt = 1e-6;
    mesh->h = 10e+12;
    mesh->volume = 1;

    mesh->temp = malloc(mesh->dim[0] * sizeof(float*));
    mesh->density = malloc(mesh->dim[0] * sizeof(float*));
    for (int i = 0; i < mesh->dim[0]; i++) {
        mesh->temp[i] = malloc(mesh->dim[1] * sizeof(float));
        mesh->density[i] = malloc(mesh->dim[1] * sizeof(float));
    }

    mesh->t_end = 1e0;

    for (int i = 0; i < mesh->dim[0]; i++) {
        for (int j = 0; j < mesh->dim[1]; j++) {
            mesh->temp[i][j] = 8e+09 + (rand() % (long)1e08);
            mesh->density[i][j] = 1e0; // Even density for debugging.
        }
    }

    flat_hydro_mesh_print(mesh);

    if (hydro_integrate_flat_mesh(mesh, options) == EXIT_FAILURE) {
        flat_hydro_mesh_destroy(&mesh);
        return EXIT_FAILURE;
    }

    flat_hydro_mesh_print(mesh);

    flat_hydro_mesh_destroy(&mesh);

    return EXIT_SUCCESS;
}

int hydro_debug(struct option_values options) {
    // linear_hydro_debug(options);
    flat_hydro_debug(options);
    return EXIT_SUCCESS;
}
