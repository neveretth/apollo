#ifndef __HIP_UTIL_H
#define __HIP_UTIL_H

#define __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>

// Return pointer to hipDeviceProp_t for best device on system. Prints info to
// stdout. Return NULL if error occurs on device creation.
struct hipDeviceProp_t* get_hip_device();

// Benchmark device, printing info to stdout. Return 0 if device fails during
// test.
int benchmark_device(struct hipDeviceProp_t* device);

#endif
