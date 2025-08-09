#ifndef __ADVANCED_OUTPUT_H
#define __ADVANCED_OUTPUT_H

#include "../../defs.h"

#include "../../types.h"

struct advout_t {
    real_t*** prev_temp;
    real_t*** prev_temp_block;
    real_t*** entropy;
    real_t*** entropy_block;
};

// Compute the timestep for entropy the given mesh.
int advout_entropy(struct advout_t* data, struct rt_hydro_mesh* mesh,
                   struct simulation_properties* sim_prop);

// Return pointer to initialized advout_t struct.
struct advout_t* advout_data_create(struct rt_hydro_mesh* mesh,
                                    struct simulation_properties sim_prop);

// Set up data for advout.
int advout_data_setup(struct advout_t* data, struct rt_hydro_mesh* mesh,
                      struct simulation_properties sim_prop);

// Free all data in advout_t structure, setting ptr to null
int advout_data_destroy(struct advout_t** data,
                        struct simulation_properties sim_prop);

#endif
