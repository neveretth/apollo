#ifndef __KERNEL_THERMONUCLEAR_H
#define __KERNEL_THERMONUCLEAR_H

#include "../../types.h"

// Integrate thermonuclear network.
int tnn_integration_kernel(real_t* P0, real_t* P1, real_t* P2, real_t* P3,
                           real_t* P4, real_t* P5, real_t* P6, real_t* Prefac,
                           real_t* Q, real_t* Rate, real_t* Flux, real_t* Fplus,
                           real_t* Fminus, real_t* FplusFac, real_t* FminusFac,
                           real_t* FplusSum, real_t* FminusSum, int* FplusMax,
                           int* FminusMax, int* MapFplus, int* MapFminus,
                           real_t* Y, int* NumReactingSpecies, int* Reactant1,
                           int* Reactant2, int* Reactant3, int number_species,
                           int number_reactions, int f_plus_total,
                           int f_minus_total, real_t t9, real_t t_max,
                           real_t dt_init);

int check_asy(real_t Fminus, real_t Y, real_t dt);

// Return the Y specified by speciesIndex updated using the forward
// Euler method.
real_t euler_update(real_t FplusSum, real_t FminusSum, real_t dt);

// Return the updated Y using the asymptotic formula.
real_t asymptotic_update(real_t Fplus, real_t Fminus, real_t Y, real_t dt);

// Return the network timestep. For now it is a placeholder
// just returning the timestep as a fixed fraction of the time.
real_t compute_timestep(real_t prevdt, real_t t, real_t tmax);

// Function to compute the effective decay constant keff for asymptotic
// approximation
// *** NOT PRESENTLY USED ***
real_t compute_keff(real_t Fminus, real_t Y);

// Perform preprocessing on thermonuclear network data.
int tnn_data_preprocess(struct tnn**** tnn, struct rt_hydro_mesh* mesh,
                        struct simulation_properties sim_prop,
                        struct option_values options);

int tnn_data_postprocess(struct tnn**** tnn, struct rt_hydro_mesh* mesh,
                        struct simulation_properties sim_prop,
                        struct option_values options);

// Update tnn kernel parameters.
int problem_parameters_update(struct problem_parameters* params,
                              struct rate_library* rates, struct tnn* network);

void reaction_mask_update(int** mask, struct rate_library* rates,
                          struct tnn* network,
                          struct problem_parameters* params, int* temp_int1,
                          int* temp_int2);

int tnn_integrate_network(struct simulation_properties sim_prop,
                          struct rate_library* rates, struct tnn* network,
                          struct problem_parameters* params,
                          struct option_values options);

int tnn_kernel_trigger(struct simulation_properties sim_prop,
                       struct rate_library* rates, struct tnn**** network,
                       struct problem_parameters* params,
                       struct option_values options);

#ifdef __MP_ROCM
#define __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>

__global__ void
integration_kernel_hip(real_t* P0, real_t* P1, real_t* P2, real_t* P3,
                       real_t* P4, real_t* P5, real_t* P6, real_t* Prefac,
                       real_t* Q, real_t* Rate, real_t* Flux, real_t* Fplus,
                       real_t* Fminus, real_t* FplusFac, real_t* FminusFac,
                       real_t* FplusSum, real_t* FminusSum, int* FplusMax,
                       int* FminusMax, int* MapFplus, int* MapFminus, real_t* Y,
                       int* NumReactingSpecies, int* Reactant1, int* Reactant2,
                       int* Reactant3, int* int_val, real_t* real_t_val);
#endif // __MP_ROCM

#endif
