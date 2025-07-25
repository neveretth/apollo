# GPU (MP)

Currently ROCm (HIP, AMD GPUs) is the only MP backed in development for 
Apollo.

## Host code
Host code segments dedicated to GPU code management are held within
preprocessor statements (`#ifdef __MP_ROCM`). These sections are omitted in 
the non-rocm (std) implementation of Apollo.

## Kernels
HIP kernels are in the kernel.hip file within a module alongside the kernel.c
serial code.
