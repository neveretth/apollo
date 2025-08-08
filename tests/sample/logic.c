#include "logic.h"
#include "../lib/lib.h"
#include "link.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int sample_kernel_data_init(void** param, int* size) {
    const int len = 100;
    size[SEQUENCE] = len * sizeof(real_t);
    param[SEQUENCE] = calloc(len, sizeof(real_t));
    // calloc initializes to 0

    size[LENGTH] = SCALAR * sizeof(int);
    param[LENGTH] = malloc(sizeof(real_t));
    *(int*)param[LENGTH] = len;

    return EXIT_SUCCESS;
}

int sample_kernel_driver(void** param, int* size) {
    sample_kernel(param[SEQUENCE], *(int*)param[LENGTH]);
    return EXIT_SUCCESS;
}

int sample_kernel_bad_driver(void** param, int* size) {
    sample_kernel_bad(param[SEQUENCE], *(int*)param[LENGTH]);
    return EXIT_SUCCESS;
}

int sample_kernel_data_check(void** param_base, int* size_base,
                             void** param_check, int* size_check) {

    const real_t tolerance = 1e-2;

    for (int i = 0; i < *(int*)param_base[LENGTH]; i++) {
        if (fabsf(((real_t*)param_base[SEQUENCE])[i] -
                  ((real_t*)param_check[SEQUENCE])[i]) > tolerance) {
            printf("[APOLLO UNIT TEST] Tolerance failure at sequence[%i]\n", i);
            printf("[APOLLO UNIT TEST] %f(read) vs %f\n",
                   ((real_t*)param_base[SEQUENCE])[i],
                   ((real_t*)param_check[SEQUENCE])[i]);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
