#ifndef __ARGS_H
#define __ARGS_H

#include "types.h"

// Returns options class for given argc, argv. This function will exit if
// certian conditions are not met.
struct option_values parse_args(int argc, char** argv);

#endif
