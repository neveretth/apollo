#ifndef __TYPES_H
#define __TYPES_H

#include "neutrino-types.h"
#include "tnn-types.h"
#include "hydro-types.h"

#include <stdio.h>
#include <H5Include.h>

// Contain values for all options when calling program.
struct option_values {
    FILE* rate_library_file;
    FILE* network_file;
    hid_t neutrino_file;
    int verbose;
    int rocm_accel;
    int rocm_debug;
    int halt;
    int thermo_debug;
    int neutrino_debug;
    int hydro_debug;
};

// Free N pointers at ptr.
int freenptr(void** ptr, int N);

// Print abundances present in network.
int print_abundances(const struct tnn* network);

// Print result data.
int print_results(const struct rate_library* rates, const struct tnn* network,
                  const struct problem_parameters* params);

// Clean up all options given passed to code, freeing pointers, clearing caches,
// etc.
int options_clean(struct option_values options);

#endif
