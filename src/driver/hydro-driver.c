#include "hydro-driver.h"
#include "kernel-driver.h"

#include <stdlib.h>

int hydro_debug(struct option_values options) {
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
