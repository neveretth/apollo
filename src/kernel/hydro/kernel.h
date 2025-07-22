#ifndef __KERNEL_HYDRO_H
#define __KERNEL_HYDRO_H

// Integrate linear hydro mesh.
int hydro_integration_kernel(float* temp, float* density, float volume, float h,
                             float dt, float t_end, int dim);

// Integrate flat hydro mesh.
int flat_hydro_integration_kernel(float** temp, float** density, float volume, float h,
                             float dt, float t_end, int* dim);

int rt_hydro_integration_kernel(float*** temp, float*** density, float volume, float h,
                             float dt, float t_end, int* dim);

#endif
