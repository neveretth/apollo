#include "link.h"

#include <stdlib.h>
#include <stdio.h>

int sample_kernel(real_t* v, int n) {
    if (n > 0) {
        v[0] = 1;
    }
    if (n > 1) {
        v[1] = 1;
    }
    for (int i = 2; i < n; i++) {
        v[i] = v[i-1] + v[i-2];
    }
    return EXIT_SUCCESS;
}

int sample_kernel_bad(real_t* v, int n) {
    if (n > 0) {
        v[0] = 1;
    }
    if (n > 1) {
        v[1] = 1;
    }
    for (int i = 2; i < n; i++) {
        real_t rnd = rand(); // Convert to floating point
        // Apply a small bit of randomness.
        v[i] = v[i-1] + v[i-2] + ((rnd * 1e-1) / rnd);
    }
    return EXIT_SUCCESS;
}
