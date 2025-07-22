#ifndef __KERNEL_NEUTRINO_H
#define __KERNEL_NEUTRINO_H

#include "../../types.h"

// Integrate neutrino network.
int neunet_integration_kernel(float** rate_in, float** rate_out, float* n_old,
                              float* ec, float* dv, float dt, float t_end,
                              float EpsA, float EpsR, float g_a, float g_b,
                              float g_c, int n_g, int halt);

// Neutrino integration function.
void compute_rates(float* F, float* k, float** R_In, float** R_Out, float* N,
                   float* dV, int n_g);

// Neutrino integration function.
void QSS2(float* Nnew, float* F0, float* Fp, float* k0, float* kp,
          float* Alpha0, float dt, float* Nold, int n_g);

// Neutrino integration function.
void QSS1(float* N_p, float* Alpha, float* F, float* k, float dt, float* Nold,
          int n_g);

// Neutrino integration function.
float compute_next_timestep(const float* E_O, const float* E_D, float dt_old,
                            int n_g);

// Perform preprocessing on neutrino network data.
int neunet_data_preprocess(struct neunet**** neunet, struct rt_hydro_mesh* mesh,
                           struct simulation_properties sim_prop);

#endif
