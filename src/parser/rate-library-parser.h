#ifndef __RATE_LIB_PARSER_H
#define __RATE_LIB_PARSER_H

#include "../types.h"

// Return pointer to rate_library struct containing rate library data. This is
// required only once, at the beginning of the entire calculation.
struct rate_library* rate_library_create(struct option_values options);

#endif
