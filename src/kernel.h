#ifndef __KERNEL_H
#define __KERNEL_H

#define __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>

// Integrate thermonuclear network.
int tnn_integration_kernel(float* P0, float* P1, float* P2, float* P3,
                           float* P4, float* P5, float* P6, float* Prefac,
                           float* Q, float* Rate, float* Flux, float* Fplus,
                           float* Fminus, float* FplusFac, float* FminusFac,
                           float* FplusSum, float* FminusSum, int* FplusMax,
                           int* FminusMax, int* MapFplus, int* MapFminus,
                           float* Y, int* NumReactingSpecies, int* Reactant1,
                           int* Reactant2, int* Reactant3, int number_species,
                           int number_reactions, int f_plus_total,
                           int f_minus_total, float t9, float t_max,
                           float dt_init, int halt);

int check_asy(float Fminus, float Y, float dt);

// Return the Y specified by speciesIndex updated using the forward
// Euler method.
float euler_update(float FplusSum, float FminusSum, float dt);

// Return the updated Y using the asymptotic formula.
float asymptotic_update(float Fplus, float Fminus, float Y, float dt);

// Return the network timestep. For now it is a placeholder
// just returning the timestep as a fixed fraction of the time.
float compute_timestep(float prevdt, float t, float tmax);

// Function to compute the effective decay constant keff for asymptotic
// approximation
// *** NOT PRESENTLY USED ***
float compute_keff(float Fminus, float Y);

// Integrate neutrino network.
int neunet_integration_kernel(float** rate_in, float** rate_out, float* n_old,
                              float* ec, float* dv, float dt, float t_end,
                              float EpsA, float EpsR, float g_a, float g_b, float g_c,
                              int n_g, int halt);

// Neutrino integration function.
void compute_rates(float* F, float* k, float** R_In, float** R_Out, float* N,
                   float* dV, int n_g);

// Neutrino integration function.
void QSS2(float* Nnew, float* F0, float* Fp, float* k0, float* kp,
          float* Alpha0, float dt, float* Nold, int n_g);

// Neutrino integration function.
void QSS1(float* N_p, float* Alpha, float* F, float* k, float dt, float* Nold, int n_g);

// Neutrino integration function.
float compute_next_timestep(const float* E_O, const float* E_D, float dt_old, int n_g);

// *****************************

// temp kernel for debugging
__global__ void vector_mult_kernel(float* A, float* B, float* C);

__global__ void integration_kernel_hip(
    float* P0, float* P1, float* P2, float* P3, float* P4, float* P5, float* P6,
    float* Prefac, float* Q, float* Rate, float* Flux, float* Fplus,
    float* Fminus, float* FplusFac, float* FminusFac, float* FplusSum,
    float* FminusSum, int* FplusMax, int* FminusMax, int* MapFplus,
    int* MapFminus, float* Y, int* NumReactingSpecies, int* Reactant1,
    int* Reactant2, int* Reactant3, int* int_val, float* float_val);

#endif
