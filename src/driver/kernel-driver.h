#ifndef __KERNEL_DRIVER_H
#define __KERNEL_DRIVER_H

#include "../types.h"

int tnn_integrate_network(struct rate_library* rates, struct tnn* network,
                          struct problem_parameters* params,
                          struct option_values options);

int neunet_integrate_network(struct neunet* network,
                             struct option_values options);

int hydro_integrate_mesh(struct rt_hydro_mesh* mesh,
                         struct option_values options);

#endif
