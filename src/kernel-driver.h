#ifndef __KERNEL_DRIVER_H
#define __KERNEL_DRIVER_H

#include "types.h"

int integrate_network(struct rate_library* rates, struct network* network,
                      struct problem_parameters* params,
                      struct option_values options);

#endif
