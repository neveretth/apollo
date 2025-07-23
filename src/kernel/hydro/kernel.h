#ifndef __KERNEL_HYDRO_H
#define __KERNEL_HYDRO_H

#include "../../defs.h"

#include "../../types.h"

// Primary hydro kernel.
int hydro_integration_kernel(real_t*** temp, real_t*** density, real_t volume, real_t h,
                             real_t dt, real_t t_end, int* dim);

// Perform preprocessing operations before the primary hydro kernel.
int hydro_data_preprocess();

// Hydro integration kernel driver.
int hydro_integrate_mesh(struct rt_hydro_mesh* mesh,
                            struct option_values options);

#endif
