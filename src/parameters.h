#ifndef __PARAMETERS_H
#define __PARAMETERS_H

#include "types.h"

// Return pointer to problem_parameters struct containing parameters for a
// simulation from the given network and rate_library.
struct problem_parameters*
problem_parameters_create(struct rate_library* rates, struct network* network,
                          struct option_values options);

// Return reaction_mask for problem_parameters_create. Alters values in the rate_library, network, and problem_parameters structs.
int** reaction_mask_create(struct rate_library* rates, struct network* network,
                           struct problem_parameters* params, struct option_values options, int* temp_int1, int* temp_int2);

#endif
