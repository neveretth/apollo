#ifndef __SAMPLE_LINK_H
#define __SAMPLE_LINK_H

// Include standardized types like real_t
#include "../../src/defs.h"

// Index order of the parameters taken by the kernel.
// These will give the correct values when subscripted.
enum SAMPLE_KERNEL_PARAMS {
  SEQUENCE = 0,
  LENGTH,
  SAMPLE_KERNEL_NUM_PARAMS // This __MUST__ be the last enum type.
};

// Normally we would include a kernel:
//
// #include "../../src/kernel/hydro/kernel.h"
//
// However for this sample we will use our own kernel, defined
// in kernel.c, a file that _only_ exists to hold this in place of 
// using a real Apollo kernel.

// Take a pointer v of length n and place the
// first n digits of the fibonacci sequence into v.
int sample_kernel(real_t* v, int n);

// This kernel operates exactly like sample_kernel, but it 
// applies randomness, causing the unit test to fail if the randomness
// falls outside the allowed tolerance.
int sample_kernel_bad(real_t* v, int n);

#endif
