#include "unified-driver.h"

#include "kernel-driver.h"
#include "../display.h"

#include <stdlib.h>

int unified_driver(struct simulation_properties sim_prop,
                   struct option_values options) {

    struct rt_hydro_mesh* mesh = malloc(sizeof(struct rt_hydro_mesh));
    mesh->dim[0] = sim_prop.resolution[0];
    mesh->dim[1] = sim_prop.resolution[1];
    mesh->dim[2] = sim_prop.resolution[2];

    // Eventually move these to TOML
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

    // INIT DATA
    for (int i = 0; i < mesh->dim[0]; i++) {
        for (int j = 0; j < mesh->dim[1]; j++) {
            for (int k = 0; k < mesh->dim[2]; k++) {
                float k_ = k;
                mesh->temp[i][j][k] =
                    8e+09 + ((16 / (k_ + 1)) * (rand() % (long)1e09));
                mesh->density[i][j][k] = 1e0; // Even density for debugging.
            }
        }
    }

    float t_end = sim_prop.t_end;
    float t = 0;
    float tres = sim_prop.output_tres;
    float t_inter = t_end / tres;

    print_float_3d(mesh->temp, mesh->dim[0], mesh->dim[1], mesh->dim[2]);
    
    while (t < t_end) {
        mesh->t_end = t_inter;
        if (sim_prop.hydro) {
            if (hydro_integrate_rt_mesh(mesh, options) == EXIT_FAILURE) {
                rt_hydro_mesh_destroy(&mesh);
                return EXIT_FAILURE;
            }
        }
        if (sim_prop.neutrino) {
            printf("not implemented\n");
            goto exit_fail;
        }
        if (sim_prop.thermo) {
            printf("not implemented\n");
            goto exit_fail;
        }
        t += t_inter;
        // printf("time: %f/%f\n", t, t_end);
    }
    printf("\n");

    print_float_3d(mesh->temp, mesh->dim[0], mesh->dim[1], mesh->dim[2]);

    rt_hydro_mesh_destroy(&mesh);
    return EXIT_SUCCESS;
exit_fail:
    rt_hydro_mesh_destroy(&mesh);
    return EXIT_FAILURE;
}
