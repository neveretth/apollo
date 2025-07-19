#include "unified-driver.h"

#include "../display.h"
#include "../numeffect.h"
#include "../parser.h"
#include "kernel-driver.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>

#define BOLTZMANN_CONSTANT 1.380658e-16

int unified_driver(struct simulation_properties sim_prop,
                   struct option_values options) {

    struct rt_hydro_mesh* mesh = malloc(sizeof(struct rt_hydro_mesh));
    struct neunet**** neutrino;
    struct tnn**** thermo;
    struct rate_library* rates;

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

    // For now hydro mesh must always be initialized, even if the compute kernel
    // is not run.
    // if (sim_prop.hydro) {
    if (1) {
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
        for (int i = 0; i < mesh->dim[0]; i++) {
            for (int j = 0; j < mesh->dim[1]; j++) {
                for (int k = 0; k < mesh->dim[2]; k++) {
                    mesh->temp[i][j][k] = base_temp;
                    mesh->density[i][j][k] = base_density;
                }
            }
        }
    }

    options.neutrino_file = sim_prop.neutrino_file;
    options.rate_library_file = sim_prop.rate_library_file;
    options.network_file = sim_prop.network_file;

    if (sim_prop.thermo) {
        rates = rate_library_create(options);
        thermo = malloc(sim_prop.resolution[0] * sizeof(struct tnn*));
        for (int i = 0; i < sim_prop.resolution[0]; i++) {
            thermo[i] = malloc(sim_prop.resolution[1] * sizeof(struct tnn*));
            for (int j = 0; j < sim_prop.resolution[1]; j++) {
                thermo[i][j] =
                    malloc(sim_prop.resolution[2] * sizeof(struct tnn*));
            }
        }
        thermo[0][0][0] = network_create(options);
        for (int i = 0; i < sim_prop.resolution[0]; i++) {
            for (int j = 0; j < sim_prop.resolution[1]; j++) {
                for (int k = 0; k < sim_prop.resolution[2]; k++) {
                    if ((i + j + k) > 0) {
                        thermo[i][j][k] = tnn_clone(thermo[0][0][0]);
                    }
                }
            }
        }
    }

    // This is proof-of-concept. It does _not_ function how a real simulation
    // should.
    if (sim_prop.neutrino) {
        neutrino = malloc(sim_prop.resolution[0] * sizeof(struct neunet*));
        for (int i = 0; i < sim_prop.resolution[0]; i++) {
            neutrino[i] =
                malloc(sim_prop.resolution[1] * sizeof(struct neunet*));
            for (int j = 0; j < sim_prop.resolution[1]; j++) {
                neutrino[i][j] =
                    malloc(sim_prop.resolution[2] * sizeof(struct neunet*));
            }
        }
        neutrino[0][0][0] = neunet_create(options);
        neutrino[0][0][0]->f->cycle_max = (int)1e9;
        neutrino[0][0][0]->f->EpsA = 1.0e-06;
        neutrino[0][0][0]->f->EpsR = 1.0e-04;
        neutrino[0][0][0]->f->tol = 0.0;
        neutrino[0][0][0]->f->g_a = 7.5e-01;
        neutrino[0][0][0]->f->g_b = 1.0e+02;
        neutrino[0][0][0]->f->g_c = sqrt(50.0);
        neutrino[0][0][0]->f->dt = mesh->dt;
        neutrino[0][0][0]->f->dt_new = neutrino[0][0][0]->f->dt;
        for (int i = 0; i < sim_prop.resolution[0]; i++) {
            for (int j = 0; j < sim_prop.resolution[1]; j++) {
                for (int k = 0; k < sim_prop.resolution[2]; k++) {
                    if ((i + j + k) > 0) {
                        neutrino[i][j][k] = neunet_clone(neutrino[0][0][0]);
                    }
                }
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

    if (sim_prop.output && sim_prop.hydro) {
        fprint_float_3d(sim_prop.hydro_out_file, mesh->temp, mesh->dim[0],
                        mesh->dim[1], mesh->dim[2]);
    }

    printf("\n\n");

    time_t kerneltime = clock();

    struct problem_parameters* params;
    while (t < t_end) {
        mesh->t_end = t_inter;
        if (sim_prop.hydro) {
            if (hydro_integrate_rt_mesh(mesh, options) == EXIT_FAILURE) {
                rt_hydro_mesh_destroy(&mesh);
                return EXIT_FAILURE;
            }
        }
        if (sim_prop.neutrino) {
            for (int i = 0; i < sim_prop.resolution[0]; i++) {
                for (int j = 0; j < sim_prop.resolution[1]; j++) {
                    for (int k = 0; k < sim_prop.resolution[2]; k++) {
                        // Hydro ==> Neutrino
                        neutrino[i][j][k]->f->kt =
                            mesh->temp[i][j][k] * BOLTZMANN_CONSTANT;
                        neutrino[i][j][k]->f->rho = mesh->density[i][j][k];
                        // TODO: regenerate certain params that change with
                        // this.
                        if (neunet_integrate_network(neutrino[i][j][k],
                                                     options) == EXIT_FAILURE) {
                            return EXIT_FAILURE;
                        }
                        // TODO: Retrieve values from neutrino.
                        // Neutrino ==> Hydro
                    }
                }
            }
        }
        if (sim_prop.thermo) {
            for (int i = 0; i < sim_prop.resolution[0]; i++) {
                for (int j = 0; j < sim_prop.resolution[1]; j++) {
                    for (int k = 0; k < sim_prop.resolution[2]; k++) {
                        // Hydro ==> Thermo
                        thermo[i][j][k]->f->rho = mesh->density[i][j][k];
                        thermo[i][j][k]->f->t9 = mesh->temp[i][j][k] / 1e9;
                        float density[3];
                        density[0] = 1.0f;
                        density[1] = thermo[i][j][k]->f->rho;
                        density[2] =
                            thermo[i][j][k]->f->rho * thermo[i][j][k]->f->rho;
                        for (int n = 0; n < rates->number_reactions; n++) {
                            rates->prefactor[n] *=
                                (density[rates->num_react_species[n] - 1]);
                        }
                        params = problem_parameters_create(
                            rates, thermo[i][j][k], options);
                        if (tnn_integrate_network(rates, thermo[i][j][k],
                                                  params,
                                                  options) == EXIT_FAILURE) {
                            return EXIT_FAILURE;
                        }
                        problem_parameters_destroy(
                            &params, thermo[i][j][k]->info->number_species);
                    }
                    // Thermo ==> Hydro
                    // TODO: Retrieve values from thermo.
                }
            }
        }
        if (sim_prop.output && sim_prop.hydro) {
            fprint_float_3d(sim_prop.hydro_out_file, mesh->temp, mesh->dim[0],
                            mesh->dim[1], mesh->dim[2]);
        }
        t += t_inter;
        printf("\x1b[1A\x1b[2K\x1b[0G  Time: [%6.2f/%6.2f]\n", t, t_end);
    }
    kerneltime = clock() - kerneltime;
    printf("\n==apollo== Simulation complete.\n");
    if (sim_prop.print_kernel_time) {
        float kerneltime_seconds = kerneltime;
        kerneltime_seconds /= CLOCKS_PER_SEC;
        printf("==apollo== Kernel clock time: %f (s)\n", kerneltime_seconds);
    }

    rt_hydro_mesh_destroy(&mesh);
    if (sim_prop.thermo) {
        rate_library_destroy(&rates);
    }
    for (int i = 0; i < sim_prop.resolution[0]; i++) {
        for (int j = 0; j < sim_prop.resolution[1]; j++) {
            for (int k = 0; k < sim_prop.resolution[2]; k++) {
                if (sim_prop.neutrino) {
                    neunet_destroy(&neutrino[i][j][k]);
                }
                if (sim_prop.thermo) {
                    network_destroy(&thermo[i][j][k]);
                }
            }
        }
    }
    return EXIT_SUCCESS;
exit_fail:
    rt_hydro_mesh_destroy(&mesh);
    return EXIT_FAILURE;
}
