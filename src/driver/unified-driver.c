#include "unified-driver.h"

#include "../display.h"
#include "../numeffect.h"
#include "kernel-driver.h"

#include <stdlib.h>

int unified_driver(struct simulation_properties sim_prop,
                   struct option_values options) {

    struct rt_hydro_mesh* mesh = malloc(sizeof(struct rt_hydro_mesh));

    // Eventually move these to TOML
    mesh->h = 10e+12;
    mesh->volume = 1;
    float t_end = sim_prop.t_end;
    mesh->dt = t_end / 100000;
    float t = 0;
    float tres = sim_prop.output_tres;
    float t_inter = t_end / tres;
    float base_temp = sim_prop.hydro_temp_base;
    float base_density = sim_prop.hydro_density_base;

    mesh->dim[0] = sim_prop.resolution[0];
    mesh->dim[1] = sim_prop.resolution[1];
    mesh->dim[2] = sim_prop.resolution[2];

    mesh->volume = mesh->volume / mesh->dim[0];
    mesh->volume = mesh->volume / mesh->dim[1];
    mesh->volume = mesh->volume / mesh->dim[2];

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
                mesh->temp[i][j][k] = base_temp;
                mesh->density[i][j][k] = base_density;
            }
        }
    }

    // APPLY SOME EFFECTS.
    if (sim_prop.hydro_temp_effect != NULL) {
        sim_prop.hydro_temp_effect(mesh->temp, mesh->dim[0], mesh->dim[1],
                                   mesh->dim[2]);
    }
    if (sim_prop.hydro_density_effect != NULL) {
        sim_prop.hydro_density_effect(mesh->density, mesh->dim[0], mesh->dim[1],
                                      mesh->dim[2]);
    }

    if (sim_prop.output) {
        fprint_float_3d(sim_prop.hydro_out_file, mesh->temp, mesh->dim[0],
                        mesh->dim[1], mesh->dim[2]);
    }

    printf("\n\n");

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
        if (sim_prop.output) {
            fprint_float_3d(sim_prop.hydro_out_file, mesh->temp, mesh->dim[0],
                            mesh->dim[1], mesh->dim[2]);
        }
        t += t_inter;
        printf("\x1b[1A\x1b[2K\x1b[0G  Time: [%6.2f/%6.2f]\n", t, t_end);
    }
    printf("\n==apollo== Simulation complete.\n");

    rt_hydro_mesh_destroy(&mesh);
    return EXIT_SUCCESS;
exit_fail:
    rt_hydro_mesh_destroy(&mesh);
    return EXIT_FAILURE;
}
