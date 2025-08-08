#ifndef __SAMPLE_LOGIC_H
#define __SAMPLE_LOGIC_H

// Driver for the sample kernel.
int sample_kernel_driver(void** param, int* size);

// Driver for the bad sample kernel.
int sample_kernel_bad_driver(void** param, int* size);

// Initialize memory for sample kernel.
int sample_kernel_data_init(void** param, int* size);

// Check param_check values for error.
int sample_kernel_data_check(void** param_base, int* size_base,
                             void** param_check, int* size_check);

#endif
