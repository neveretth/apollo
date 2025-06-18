#ifndef __KERNEL_H
#define __KERNEL_H

int integration_kernel(float* P0, float* P1, float* P2, float* P3, float* P4,
                       float* P5, float* P6, float* Prefac, float* Q,
                       float* Rate, float* Flux, float* Fplus, float* Fminus,
                       float* FplusFac, float* FminusFac, float* FplusSum,
                       float* FminusSum, int* FplusMax, int* FminusMax,
                       int* MapFplus, int* MapFminus, float* Y,
                       int* NumReactingSpecies, int* Reactant1, int* Reactant2,
                       int* Reactant3, int number_species, int number_reactions, int f_plus_total, int f_minus_total, float t9, float t_max, float dt_init, int halt);

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

#endif
