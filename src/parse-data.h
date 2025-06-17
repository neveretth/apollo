#ifndef __PARSE_DATA_H
#define __PARSE_DATA_H

#include "types.h"

// Return pointer to rate_library struct containing rate library data. This is
// required only once, at the beginning of the entire calculation.
struct rate_library* rate_library_create(struct option_values options);

// Return pointer to network struct containing associated partition functions.
// This is required only once at the very beginning of the calculation.
struct network* network_create(struct option_values options);

#endif
