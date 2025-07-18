#include "hydro-driver.h"
#include "../display.h"
#include "kernel-driver.h"

#include <stdlib.h>
#include <math.h>

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

// Debug code for generating a graph (heatmap) of hydro evolution.
static int flat_hydro_graph_debug(struct option_values options) {
    struct flat_hydro_mesh* mesh = malloc(sizeof(struct flat_hydro_mesh));
    mesh->dim[0] = 256;
    mesh->dim[1] = 256;

    mesh->dt = 1e-8;
    mesh->h = 10e+12;
    mesh->volume = 1;

    mesh->temp = malloc(mesh->dim[0] * sizeof(float*));
    mesh->density = malloc(mesh->dim[0] * sizeof(float*));
    for (int i = 0; i < mesh->dim[0]; i++) {
        mesh->temp[i] = malloc(mesh->dim[1] * sizeof(float));
        mesh->density[i] = malloc(mesh->dim[1] * sizeof(float));
    }

    for (int i = 0; i < mesh->dim[0]; i++) {
        for (int j = 0; j < mesh->dim[1]; j++) {
            mesh->temp[i][j] = 8e+09 + (rand() % (long)1e09);
            mesh->density[i][j] = 1e0; // Even density for debugging.
        }
    }

    float t_end = 1e-4;
    float t = 0;
    float t_inter = t_end / 100;

    while (t < t_end) {
        // Recall that integration kernels always assume the starting time is 0.
        mesh->t_end = t_inter;
        print_float_2d(mesh->temp, mesh->dim[0], mesh->dim[1]);
        if (hydro_integrate_flat_mesh(mesh, options) == EXIT_FAILURE) {
            flat_hydro_mesh_destroy(&mesh);
            return EXIT_FAILURE;
        }
        t += t_inter;
    }

    print_float_2d(mesh->temp, mesh->dim[0], mesh->dim[1]);

    flat_hydro_mesh_destroy(&mesh);

    return EXIT_SUCCESS;
}

static int rt_hydro_graph_debug(struct option_values options) {
    struct rt_hydro_mesh* mesh = malloc(sizeof(struct rt_hydro_mesh));
    mesh->dim[0] = 8;
    mesh->dim[1] = 8;
    mesh->dim[2] = 8;

    mesh->dt = 1e-4;
    mesh->h = 10e+10;
    mesh->volume = 1;

    mesh->temp = malloc(mesh->dim[0] * sizeof(float*));
    mesh->density = malloc(mesh->dim[0] * sizeof(float*));
    for (int i = 0; i < mesh->dim[0]; i++) {
        mesh->temp[i] = malloc(mesh->dim[1] * sizeof(float*));
        mesh->density[i] = malloc(mesh->dim[1] * sizeof(float*));
        for (int j = 0; j < mesh->dim[1]; j++) {
            mesh->temp[i][j] = malloc(mesh->dim[2] * sizeof(float));
            mesh->density[i][j] = malloc(mesh->dim[2] * sizeof(float));
        }
    }

    for (int i = 0; i < mesh->dim[0]; i++) {
        for (int j = 0; j < mesh->dim[1]; j++) {
            for (int k = 0; k < mesh->dim[2]; k++) {
                float k_ = k;
                mesh->temp[i][j][k] = 8e+09 + ((16 / (k_ + 1)) * (rand() % (long)1e09));
                mesh->density[i][j][k] = 1e0; // Even density for debugging.
            }
        }
    }

    float t_end = 1e00;
    float t = 0;
    float t_inter = t_end / 800;

    while (t < t_end) {
        // Recall that integration kernels always assume the starting time is 0.
        mesh->t_end = t_inter;
        // print_float_3d(mesh->temp, mesh->dim[0], mesh->dim[1], mesh->dim[2]);
        // printf("time: %f/%f\n", t, t_end);
        if (hydro_integrate_rt_mesh(mesh, options) == EXIT_FAILURE) {
            rt_hydro_mesh_destroy(&mesh);
            return EXIT_FAILURE;
        }
        t += t_inter;
    }
    printf("\n");

    print_float_3d(mesh->temp, mesh->dim[0], mesh->dim[1], mesh->dim[2]);

    rt_hydro_mesh_destroy(&mesh);

    return EXIT_SUCCESS;
}

int hydro_debug(struct option_values options) {
    // linear_hydro_debug(options);
    // flat_hydro_debug(options);
    // flat_hydro_graph_debug(options);
    rt_hydro_graph_debug(options);
    return EXIT_SUCCESS;
}
