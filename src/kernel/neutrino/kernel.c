#include "kernel.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef __MP_ROCM
#include "../../rocm/hip-util.h"
#endif

enum NEUTRINO_REAL_VAL {
    NEUREALVAL_DT = 0,
    NEUREALVAL_TEMP,
    NEUREALVAL_NUM_ARGS,
};

// Migrate this to real_t_val and int_val (as the GPU kernel does).
int neunet_integration_kernel(real_t* real_val, int* int_val,
                              real_t** rate_in, real_t** rate_out,
                              real_t* n_old, real_t* ec, real_t* dv, real_t dt,
                              real_t t_end, real_t EpsA, real_t EpsR,
                              real_t g_a, real_t g_b, real_t g_c, int n_g) {

    // Initialize n_eq based on mu and kt for the current model
    // real_t* n_eq = malloc(sizeof(real_t) * n_g);
    // for (int i = 0; i < n_g; i++) {
    //     n_eq[i] =
    //         1.0 /
    //         (exp((ec[i] - mu) / kt) +
    //          1.0);
    // }

    dt = real_val[NEUREALVAL_DT];
    real_t kt = real_val[NEUREALVAL_TEMP];
    
    for (int i = 0; i < n_g; i++) {
        for (int j = 0; j < n_g; j++) {
            if (j < i) {
                rate_in[i][j] = rate_in[j][i] * exp((ec[i]) / kt);
            }
        }
    }

    for (int i = 0; i < n_g; i++) {
        for (int j = 0; j < n_g; j++) {
            rate_out[i][j] = rate_in[j][i];
        }
    }

    real_t* n_0 = malloc(sizeof(real_t) * n_g);

    for (int i = 0; i < n_g; i++) {
        n_0[i] = g_a * exp(-0.5 * pow((ec[i] - g_b) / g_c, 2));
    }

    for (int i = 0; i < n_g; i++) {
        n_old[i] = n_0[i];
    }

    int done = 0;
    int restep = 0;
    int cycle = 0;
    int true_cycle = 0;
    real_t dt_new = dt;
    real_t t = 0;

    real_t* n_new = malloc(sizeof(real_t) * n_g);
    real_t* n_1 = malloc(sizeof(real_t) * n_g);
    real_t* n_2 = malloc(sizeof(real_t) * n_g);
    real_t* k0 = malloc(sizeof(real_t) * n_g);
    real_t* k1 = malloc(sizeof(real_t) * n_g);
    real_t* f0 = malloc(sizeof(real_t) * n_g);
    real_t* f1 = malloc(sizeof(real_t) * n_g);
    real_t* fp = malloc(sizeof(real_t) * n_g);
    real_t* kp = malloc(sizeof(real_t) * n_g);
    real_t* np = malloc(sizeof(real_t) * n_g);
    real_t* alpha_0 = malloc(sizeof(real_t) * n_g);
    real_t* error_desired = malloc(sizeof(real_t) * n_g);
    real_t* error_observed = malloc(sizeof(real_t) * n_g);

    while (done == 0) {
        true_cycle++;

        if (!restep) {
            cycle++;
        }

        memset(n_new, 0, n_g);

        dt = dt_new;

        compute_rates(f0, k0, rate_in, rate_out, n_old, dv, n_g);

        QSS1(np, alpha_0, f0, k0, 0.5 * dt, n_old, n_g);

        compute_rates(fp, kp, rate_in, rate_out, np, dv, n_g);

        QSS2(n_1, f0, fp, k0, kp, alpha_0, 0.5 * dt, n_old, n_g);

        compute_rates(f1, k1, rate_in, rate_out, n_old, dv, n_g);

        QSS1(np, alpha_0, f1, k1, 0.5 * dt, n_1, n_g);

        compute_rates(fp, kp, rate_in, rate_out, np, dv, n_g);

        QSS2(n_1, f0, fp, k0, kp, alpha_0, 0.5 * dt, n_1, n_g);

        QSS1(np, alpha_0, f0, k0, dt, n_old, n_g);

        compute_rates(fp, kp, rate_in, rate_out, np, dv, n_g);

        QSS2(n_2, f0, fp, k0, kp, alpha_0, dt, n_old, n_g);

        for (int i = 0; i < n_g; i++) {
            n_new[i] = n_1[i];
        }

        for (int i = 0; i < n_g; i++) {
            error_observed[i] = fabs(n_1[i] - n_2[i]);
            error_desired[i] = EpsA + EpsR * fabs(n_1[i]);
        }

        dt_new = compute_next_timestep(error_observed, error_desired, dt, n_g);
        if (dt_new > (t_end - t)) {
            dt_new = t_end - t;
        }

        for (int i = 0; i < n_g; i++) {
            n_old[i] = n_new[i];
        }

        t += dt;

        if ((t - t_end) >= 0) {
            done = 1;
            break;
        }
    }

    // Corresponds to the opposite expression above.
    real_val[NEUREALVAL_DT] = dt;

exit:
    free(n_0);
    free(n_new);
    free(n_1);
    free(n_2);
    free(k0);
    free(k1);
    free(f0);
    free(f1);
    free(fp);
    free(kp);
    free(np);
    free(alpha_0);
    free(error_desired);
    free(error_observed);

    return EXIT_SUCCESS;
}

void QSS1(real_t* N_p, real_t* Alpha, real_t* F, real_t* k, real_t dt,
          real_t* n_old, int n_g) {
    for (int i = 0; i < n_g; i++) {
        N_p[i] = n_old[i];
        real_t r;

        r = 1.0 / (k[i] * dt);

        Alpha[i] = ((160 * pow(r, (real_t)3)) + (60 * pow(r, (real_t)2)) +
                    (11 * r) + 1) /
                   ((360 * pow(r, (real_t)3)) + (60 * pow(r, (real_t)2)) +
                    (12 * r) + 1);

        N_p[i] +=
            dt * (F[i] - (k[i] * n_old[i])) / (1.0 + (Alpha[i] * k[i] * dt));
    }
}

void QSS2(real_t* Nnew, real_t* f0, real_t* fp, real_t* k0, real_t* kp,
          real_t* alpha_0, real_t dt, real_t* n_old, int n_g) {
    for (int i = 0; i < n_g; i++) {
        Nnew[i] = n_old[i];
        // real_t Ft;
        real_t kBAR;
        real_t rBAR;
        real_t AlphaBAR;

        kBAR = 0.5 * (kp[i] + k0[i]);
        rBAR = 1.0 / (kBAR * dt);
        AlphaBAR =
            ((160 * pow(rBAR, 3)) + (60 * pow(rBAR, 2)) + (11 * rBAR) + 1) /
            ((360 * pow(rBAR, 3)) + (60 * pow(rBAR, 2)) + (12 * rBAR) + 1);
        // Ft = AlphaBAR * fp[i] + (1.0 - AlphaBAR) * f0[i];
        Nnew[i] +=
            dt * (fp[i] - kBAR * n_old[i]) / (1.0 + (alpha_0[i] * kBAR * dt));
    }
}

void compute_rates(real_t* F, real_t* k, real_t** R_In, real_t** R_Out,
                   real_t* N, real_t* dv, int n_g) {

    for (int i = 0; i < n_g; i++) {
        k[i] = 0.0;
        F[i] = 0.0;
        for (int j = 0; j < n_g; j++) {
            F[i] += R_In[i][j] * N[j] * dv[j];
            k[i] +=
                R_In[i][j] * N[j] * dv[j] + R_Out[i][j] * (1.0 - N[j]) * dv[j];
        }
    }
}

// "Temporary" for MAX/MIN
#include <sys/param.h>

real_t compute_next_timestep(const real_t* E_O, const real_t* E_D,
                             real_t dt_old, int n_g) {

    real_t E_R = 0.0;
    for (int i = 0; i < n_g; i++) {
        E_R = MAX(E_R, E_O[i] / E_D[i]);
    }

    real_t dt_new;
    if (E_R > dt_tol_high) {
        dt_new = dt_old * 0.9 * pow(E_R, -1.0);
    } else if (E_R < dt_tol_low) {
        dt_new = dt_old * 0.9 * pow(E_R, -0.5);
    } else {
        dt_new = dt_old;
    }

    dt_new = MAX(MIN(dt_new, 2.0 * dt_old), 0.5 * dt_old);

    return dt_new;
}

enum neunet_kernel_args {
    REAL_T_VAL = 0,
    INT_VAL,
    R_IN,
    R_OUT,
    EC,
    DV,
    N_EQ,
    NOLD,
    NUM_ARGS // THIS MUST BE THE LAST ONE
};

struct neunet_gpu_status {
    bool initialized;
    void***** args; // scary =ω=
    void* stream;
};

static struct neunet_gpu_status gpu_stat = {false};

#ifdef __MP_ROCM
static int device_init_args(struct neunet**** neunet,
                            struct simulation_properties sim_prop,
                            struct option_values opts) {

    // Lord please forgive me.
    gpu_stat.args = malloc(sim_prop.resolution[0] * sizeof(void*));
    for (int i = 0; i < sim_prop.resolution[0]; i++) {
        gpu_stat.args[i] = malloc(sim_prop.resolution[0] * sizeof(void*));
        for (int j = 0; j < sim_prop.resolution[1]; j++) {
            gpu_stat.args[i][j] =
                malloc(sim_prop.resolution[1] * sizeof(void*));
            for (int k = 0; k < sim_prop.resolution[2]; k++) {
                gpu_stat.args[i][j][k] = malloc(NUM_ARGS * sizeof(void*));
                for (int z = 0; z < NUM_ARGS; z++) {
                    gpu_stat.args[i][j][k][z] = malloc(sizeof(void*));
                }
            }
        }
    }

    // Create unique trackers.
    for (int i = 0; i < sim_prop.resolution[0]; i++) {
        for (int j = 0; j < sim_prop.resolution[1]; j++) {
            for (int k = 0; k < sim_prop.resolution[2]; k++) {
                if (devbuf_create(gpu_stat.args[i][j][k][EC],
                                  neunet[i][j][k]->info->num_groups *
                                      sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_create(gpu_stat.args[i][j][k][DV],
                                  neunet[i][j][k]->info->num_groups *
                                      sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_create(gpu_stat.args[i][j][k][NOLD],
                                  neunet[i][j][k]->info->num_groups *
                                      sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_create(gpu_stat.args[i][j][k][N_EQ],
                                  neunet[i][j][k]->info->num_groups *
                                      sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_create(gpu_stat.args[i][j][k][REAL_T_VAL],
                                  16 * sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_create(gpu_stat.args[i][j][k][INT_VAL],
                                  16 * sizeof(int)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                // These are flattened 2D arrays
                if (devbuf_create(gpu_stat.args[i][j][k][R_IN],
                                  (pow(neunet[i][j][k]->info->num_groups, 2)) *
                                      sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_create(gpu_stat.args[i][j][k][R_OUT],
                                  (pow(neunet[i][j][k]->info->num_groups, 2)) *
                                      sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
            }
        }
    }

    // Finally, send initial data to GPU
    for (int i = 0; i < sim_prop.resolution[0]; i++) {
        for (int j = 0; j < sim_prop.resolution[1]; j++) {
            for (int k = 0; k < sim_prop.resolution[2]; k++) {
                if (devbuf_write(gpu_stat.args[i][j][k][EC],
                                 neunet[i][j][k]->fptr->ec,
                                 neunet[i][j][k]->info->num_groups *
                                     sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_write(gpu_stat.args[i][j][k][DV],
                                 neunet[i][j][k]->fptr->dv,
                                 neunet[i][j][k]->info->num_groups *
                                     sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_write(gpu_stat.args[i][j][k][NOLD],
                                 neunet[i][j][k]->fptr->n_old,
                                 neunet[i][j][k]->info->num_groups *
                                     sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_write(gpu_stat.args[i][j][k][N_EQ],
                                 neunet[i][j][k]->fptr->n_eq,
                                 neunet[i][j][k]->info->num_groups *
                                     sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_write(gpu_stat.args[i][j][k][N_EQ],
                                 neunet[i][j][k]->fptr->n_eq,
                                 neunet[i][j][k]->info->num_groups *
                                     sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_write_flatten(gpu_stat.args[i][j][k][R_IN],
                                         (void*)neunet[i][j][k]->info->rate_in,
                                         neunet[i][j][k]->info->num_groups,
                                         sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
                if (devbuf_write_flatten(gpu_stat.args[i][j][k][R_OUT],
                                         (void*)neunet[i][j][k]->info->rate_out,
                                         neunet[i][j][k]->info->num_groups,
                                         sizeof(real_t)) != EXIT_SUCCESS) {
                    return EXIT_FAILURE;
                }
            }
        }
    }

    return EXIT_SUCCESS;
}
#endif

// This really highlights the weakness of storing the neutrino network data like
// this. A rework could be extremely beneficial for parallelization.
int neunet_data_preprocess(struct neunet**** neunet, struct rt_hydro_mesh* mesh,
                           struct simulation_properties sim_prop,
                           struct option_values opts) {
#ifdef __MP_ROCM
    if (!gpu_stat.initialized) {
        if (device_init_args(neunet, sim_prop, opts) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
        gpu_stat.initialized = true;
    }
#endif
    for (int i = 0; i < sim_prop.resolution[XDIM]; i++) {
        for (int j = 0; j < sim_prop.resolution[YDIM]; j++) {
            for (int k = 0; k < sim_prop.resolution[ZDIM]; k++) {
                // TODO: we don't know the exact units.
                // MeV?
                neunet[i][j][k]->f->kt =
                    mesh->temp[k][j][i] * 8.621738e-5 * BOLTZMANN_CONSTANT * 1e12;
                neunet[i][j][k]->f->rho = mesh->density[k][j][i];
                neunet[i][j][k]->f->t_end = mesh->dt;
            }
        }
    }

#ifdef __MP_ROCM
    // This is not ideal with how streams _should_ work.
    // This part updates params...

    real_t* tmpreal = malloc(16 * sizeof(real_t));
    int* tmpint = malloc(16 * sizeof(int));
    for (int i = 0; i < sim_prop.resolution[0]; i++) {
        for (int j = 0; j < sim_prop.resolution[1]; j++) {
            for (int k = 0; k < sim_prop.resolution[2]; k++) {
                tmpreal[0] = neunet[i][j][k]->f->t; // Should just be 0
                tmpreal[1] = neunet[i][j][k]->f->dt;
                tmpreal[2] = neunet[i][j][k]->f->g_a;
                tmpreal[3] = neunet[i][j][k]->f->g_b;
                tmpreal[4] = neunet[i][j][k]->f->g_c;
                tmpreal[5] = neunet[i][j][k]->f->tol;
                tmpreal[6] = neunet[i][j][k]->f->EpsA;
                tmpreal[7] = neunet[i][j][k]->f->EpsR;
                tmpreal[8] = neunet[i][j][k]->f->tol;
                tmpreal[9] = neunet[i][j][k]->f->err;
                tmpreal[10] = neunet[i][j][k]->f->t_end;
                tmpint[0] = neunet[i][j][k]->info->num_groups;
                devbuf_write(gpu_stat.args[i][j][k][REAL_T_VAL], tmpreal,
                             16 * sizeof(real_t));
                devbuf_write(gpu_stat.args[i][j][k][INT_VAL], tmpint,
                             16 * sizeof(int));
            }
        }
    }
    free(tmpreal);
    free(tmpint);
#endif
    return EXIT_SUCCESS;
}

int neunet_integrate_network(struct simulation_properties sim_prop,
                             struct neunet* network,
                             struct option_values options, int i, int j,
                             int k) {

#ifdef __MP_ROCM
    struct dim3 blockdim = {network->info->num_groups, 1, 1};
    struct dim3 griddim = {1, 1, 1};
    if (options.rocm_accel) {
        // Pass 0 for sharedmem because I don't know what it does?
        if (hipLaunchKernel(qss_hip_neunet_kernel, griddim, blockdim,
                            gpu_stat.args[i][j][k], 0,
                            hipStreamDefault) != hipSuccess) {
            return EXIT_FAILURE;
        }
    } else {
        if (neunet_integration_kernel(
                network->info->rate_in, network->info->rate_out,
                network->fptr->n_old, network->fptr->ec, network->fptr->dv,
                network->f->dt, network->f->t_end, network->f->EpsA,
                network->f->EpsR, network->f->g_a, network->f->g_b,
                network->f->g_c, network->info->num_groups) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
    }
#else
    // What is tmp?
    // It represents the pointer that may be introduced at network->real
    // It is currently being used in the kernel to track values persistent
    // across multiple invocations.
    real_t* tmp = malloc(NEUREALVAL_NUM_ARGS * sizeof(real_t));
    tmp[NEUREALVAL_TEMP] = network->f->kt;
    tmp[NEUREALVAL_DT] = network->f->dt;
    if (neunet_integration_kernel(
            tmp, NULL, network->info->rate_in, network->info->rate_out,
            network->fptr->n_old, network->fptr->ec, network->fptr->dv,
            network->f->dt, network->f->t_end, network->f->EpsA,
            network->f->EpsR, network->f->g_a, network->f->g_b, network->f->g_c,
            network->info->num_groups) == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }
    network->f->dt = tmp[NEUREALVAL_DT];
    free(tmp);
#endif
    return EXIT_SUCCESS;
}

int neunet_kernel_trigger(struct simulation_properties sim_prop,
                          struct neunet**** network,
                          struct option_values options) {
#ifdef __MP_ROCM
    // Something something stream....?
#endif
    for (int i = 0; i < sim_prop.resolution[0]; i++) {
        for (int j = 0; j < sim_prop.resolution[1]; j++) {
            for (int k = 0; k < sim_prop.resolution[2]; k++) {
                if (neunet_integrate_network(sim_prop, network[i][j][k],
                                             options, i, j,
                                             k) == EXIT_FAILURE) {
                    return EXIT_FAILURE;
                }
            }
        }
    }
#ifdef __MP_ROCM
    // I'm not sure how this works...
    hipStreamSynchronize(hipStreamDefault);
#endif
    return EXIT_SUCCESS;
}

int neunet_data_postprocess(struct neunet**** neunet,
                            struct rt_hydro_mesh* mesh,
                            struct simulation_properties sim_prop,
                            struct option_values opts) {
    // (READ FROM GPU AND) SEND CHANGE TO HYDRO
    return EXIT_SUCCESS;
}
