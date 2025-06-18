#include <hip/hip_runtime.h>

extern "C" __global__ void vector_mult_kernel(float* A, float* B, float* C) {
    int tid = threadIdx.x;
    C[tid] = A[tid] * B[tid];
}
