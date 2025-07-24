#ifndef __HIP_UTIL_H
#define __HIP_UTIL_H

#define __HIP_PLATFORM_AMD__
#include <hip/hip_runtime.h>

#include "../defs.h"

// Return pointer to hipDeviceProp_t for best device on system. Prints info to
// stdout. Return NULL if error occurs on device creation.
struct hipDeviceProp_t* get_hip_device();

// Benchmark device, printing info to stdout. Return 0 if device fails during
// test.
int benchmark_device(struct hipDeviceProp_t* device);

int devbuf_create(void** devptr, int size);

int devbuf_write(void** devptr, void* hostptr, int size);

// NOTE: len should be the number of rows/columns (they MUST BE EQUAL).
// size should be the size of the datatype
int devbuf_write_flatten(void** devptr, void** hostptr, int len, int size);

int devbuf_read(void* hostptr, void** devptr, int size);

#endif
