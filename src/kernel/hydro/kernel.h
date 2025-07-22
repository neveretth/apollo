#ifndef __KERNEL_HYDRO_H
#define __KERNEL_HYDRO_H

// Primary hydro kernel.
int hydro_integration_kernel(float*** temp, float*** density, float volume, float h,
                             float dt, float t_end, int* dim);

// Perform preprocessing operations before the primary hydro kernel.
int hydro_data_preprocess();

#endif
