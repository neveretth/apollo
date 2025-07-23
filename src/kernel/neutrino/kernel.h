#ifndef __KERNEL_NEUTRINO_H
#define __KERNEL_NEUTRINO_H

#include "../../types.h"

// Integrate neutrino network.
int neunet_integration_kernel(real_t** rate_in, real_t** rate_out, real_t* n_old,
                              real_t* ec, real_t* dv, real_t dt, real_t t_end,
                              real_t EpsA, real_t EpsR, real_t g_a, real_t g_b,
                              real_t g_c, int n_g);

// Neutrino integration function.
void compute_rates(real_t* F, real_t* k, real_t** R_In, real_t** R_Out, real_t* N,
                   real_t* dV, int n_g);

// Neutrino integration function.
void QSS2(real_t* Nnew, real_t* F0, real_t* Fp, real_t* k0, real_t* kp,
          real_t* Alpha0, real_t dt, real_t* Nold, int n_g);

// Neutrino integration function.
void QSS1(real_t* N_p, real_t* Alpha, real_t* F, real_t* k, real_t dt, real_t* Nold,
          int n_g);

// Neutrino integration function.
real_t compute_next_timestep(const real_t* E_O, const real_t* E_D, real_t dt_old,
                            int n_g);

// Perform preprocessing on neutrino network data.
int neunet_data_preprocess(struct neunet**** neunet, struct rt_hydro_mesh* mesh,
                           struct simulation_properties sim_prop);

int neunet_integrate_network(struct neunet* network,
                             struct option_values options);

#endif
