#include "unified-driver.h"

#include "../display.h"
#include "../numeffect.h"
#include "../parser.h"
#include "../parser/tnn-parameters.h"

#include "../kernel/hydro/kernel.h"
#include "../kernel/neutrino/kernel.h"
#include "../kernel/thermonuclear/kernel.h"

#ifdef __MP_ROCM
#include "../rocm/hip-util.h"
#endif

#include <math.h>
#include <stdlib.h>
#include <time.h>

int unified_driver(struct simulation_properties sim_prop,
                   struct option_values options) {

    struct rt_hydro_mesh* mesh = malloc(sizeof(struct rt_hydro_mesh));
    struct neunet**** neutrino;
    struct tnn**** thermo;
    struct rate_library* rates;

    // Eventually move these to TOML
    mesh->h = sim_prop.h;
    mesh->volume = sim_prop.volume;
    real_t t_end = sim_prop.t_end;
    real_t dt = sim_prop.dt_init;
    real_t t = 0;
    real_t tres = sim_prop.output_tres;
    real_t t_inter_lvl = t_end / tres;
    real_t t_inter = 0;
    real_t base_temp = sim_prop.hydro_temp_base;
    real_t base_density = sim_prop.hydro_density_base;

    // For now hydro mesh must always be initialized, even if the compute kernel
    // is not run.
    // if (sim_prop.hydro) {
    if (1) {
        mesh->dim[0] = sim_prop.resolution[0];
        mesh->dim[1] = sim_prop.resolution[1];
        mesh->dim[2] = sim_prop.resolution[2];
        mesh->dt = dt;

        mesh->volume = mesh->volume / mesh->dim[0];
        mesh->volume = mesh->volume / mesh->dim[1];
        mesh->volume = mesh->volume / mesh->dim[2];

        mesh->temp = malloc(mesh->dim[0] * sizeof(real_t*));
        mesh->density = malloc(mesh->dim[0] * sizeof(real_t*));
        for (int i = 0; i < mesh->dim[0]; i++) {
            mesh->temp[i] = malloc(mesh->dim[1] * sizeof(real_t*));
            mesh->density[i] = malloc(mesh->dim[1] * sizeof(real_t*));
            for (int j = 0; j < mesh->dim[1]; j++) {
                mesh->temp[i][j] = malloc(mesh->dim[2] * sizeof(real_t));
                mesh->density[i][j] = malloc(mesh->dim[2] * sizeof(real_t));
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

    struct problem_parameters* params;
    if (sim_prop.thermo) {
        rates = rate_library_create(sim_prop);
        thermo = malloc(sim_prop.resolution[0] * sizeof(struct tnn*));
        for (int i = 0; i < sim_prop.resolution[0]; i++) {
            thermo[i] = malloc(sim_prop.resolution[1] * sizeof(struct tnn*));
            for (int j = 0; j < sim_prop.resolution[1]; j++) {
                thermo[i][j] =
                    malloc(sim_prop.resolution[2] * sizeof(struct tnn*));
            }
        }
        thermo[0][0][0] = network_create(sim_prop);
        thermo[0][0][0]->f->dt_init = 1e-17;
        for (int i = 0; i < sim_prop.resolution[0]; i++) {
            for (int j = 0; j < sim_prop.resolution[1]; j++) {
                for (int k = 0; k < sim_prop.resolution[2]; k++) {
                    if ((i + j + k) > 0) {
                        thermo[i][j][k] = tnn_clone(thermo[0][0][0]);
                    }
                }
            }
        }
        params = problem_parameters_create(rates, thermo[0][0][0], options);
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
        neutrino[0][0][0] = neunet_create(sim_prop);
        neutrino[0][0][0]->f->cycle_max = (int)1e9;
        neutrino[0][0][0]->f->EpsA = 1.0e-06;
        neutrino[0][0][0]->f->EpsR = 1.0e-04;
        neutrino[0][0][0]->f->tol = 0.0;
        neutrino[0][0][0]->f->g_a = 7.5e-01;
        neutrino[0][0][0]->f->g_b = 1.0e+02;
        neutrino[0][0][0]->f->g_c = sqrt(50.0);
        neutrino[0][0][0]->f->dt = 1e-15;
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

    if (sim_prop.output) {
        fprint_real_t_3d(sim_prop.hydro_out_file, mesh->temp, mesh->dim[0],
                         mesh->dim[1], mesh->dim[2]);
    }

    // INIT ROCM IF APPROPRIATE
#ifdef __MP_ROCM
    if (options.rocm_accel) {
        options.rocm_device = get_hip_device();
    }
#endif

    printf("\n\n");

    time_t kerneltime = clock();

    bool fail = false;

    while (t < t_end) {
        t_inter += t_inter_lvl;
        while (t < t_inter) {
            // If the mesh adjusts the timestep.
            dt = mesh->dt;
            mesh->t_end = 0; // Single step
            if (sim_prop.hydro) {
                if (hydro_data_preprocess() == EXIT_FAILURE) {
                    fail = true;
                    goto exit;
                }
                // Not a trigger for now, probably will be at some point.
                if (hydro_integrate_mesh(sim_prop, mesh, options) ==
                    EXIT_FAILURE) {
                    fail = true;
                    goto exit;
                }
                if (hydro_data_postprocess() == EXIT_FAILURE) {
                    fail = true;
                    goto exit;
                }
            }
            if (sim_prop.neutrino) {
                if (neunet_data_preprocess(neutrino, mesh, sim_prop, options) ==
                    EXIT_FAILURE) {
                    fail = true;
                    goto exit;
                }
                if (neunet_kernel_trigger(sim_prop, neutrino, options) ==
                    EXIT_FAILURE) {
                    fail = true;
                    goto exit;
                }
                if (neunet_data_postprocess(neutrino, mesh, sim_prop,
                                            options) == EXIT_FAILURE) {
                    fail = true;
                    goto exit;
                }
            }
            if (sim_prop.thermo) {
                if (tnn_data_preprocess(thermo, mesh, sim_prop, options) ==
                    EXIT_FAILURE) {
                    fail = true;
                    goto exit;
                }
                if (tnn_kernel_trigger(sim_prop, rates, thermo, params,
                                       options) == EXIT_FAILURE) {
                    fail = true;
                    goto exit;
                }
                if (tnn_data_postprocess(thermo, mesh, sim_prop, options) ==
                    EXIT_FAILURE) {
                    fail = true;
                    goto exit;
                }
            }

            t += dt;
        }
        
        if (sim_prop.output) {
            fprint_real_t_3d(sim_prop.hydro_out_file, mesh->temp,
                             mesh->dim[0], mesh->dim[1], mesh->dim[2]);
        }
        
        printf("\x1b[1A\x1b[2K\x1b[0G  Time: [%6.2f/%6.2f]\n", t, t_end);
        // printf("%.10f\n", mesh->temp[0][0][0]);
    }

    kerneltime = clock() - kerneltime;
    printf("\n==apollo== Simulation complete.\n");
    if (sim_prop.print_kernel_time) {
        real_t kerneltime_seconds = kerneltime;
        kerneltime_seconds /= CLOCKS_PER_SEC;
        printf("==apollo== Kernel clock time: %f (s)\n", kerneltime_seconds);
    }

exit:
    rt_hydro_mesh_destroy(&mesh);
    if (sim_prop.thermo) {
        rate_library_destroy(&rates);
    }
    // Proper destruction should occur within the kernel modules...
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
    if (fail) {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
