#ifndef __PARSE_DATA_H
#define __PARSE_DATA_H

#include "../types.h"

// Return pointer to network struct containing the component network structures.
// This is required only once at the very beginning of the calculation.
struct tnn* network_create(struct simulation_properties sim_prop);

#endif
