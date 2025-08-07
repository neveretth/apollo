#ifndef __TYPES_H
#define __TYPES_H

#include "defs.h"

#include "toml/tomlc17.h"
#include "types/hydro-types.h"
#include "types/neutrino-types.h"
#include "types/tnn-types.h"

#include <H5Include.h>
#include <stdio.h>

#include <stdbool.h>

// Contain values for all options when calling program.
// Includes values from config.toml
struct option_values {
    char* config_file;
    char* simulation_file;
    char* root_dir;
    int verbose;
    bool rocm_accel;
    void* rocm_device; // Ignored if not __MP_ROCM (type is hipDeviceProp_t*)
};

enum TIMESCALE {
    TIMESCALE_LINEAR,
    TIMESCALE_LOGSKEW,
    TIMESCALE_LOG2,
    TIMESCALE_LOG10,
};

struct simulation_properties {
    int resolution[3];
    bool hydro;
    bool thermo;
    bool neutrino;
    FILE* rate_library_file;
    FILE* network_file;
    hid_t neutrino_file;
    bool output;
    bool print_kernel_time;
    real_t t_end;
    int output_tres;
    int timescale;
    FILE* temp_out_file;
    FILE* density_out_file;
    FILE* entropy_out_file;
    int (*hydro_temp_effect)(real_t*** data, int i_, int j_, int k_);
    int (*hydro_density_effect)(real_t*** data, int i_, int j_, int k_);
    real_t hydro_temp_base;
    real_t hydro_density_base;
    real_t dt_init;
    real_t volume;
    real_t h;
};

// Free N pointers at ptr.
int freenptr(void** ptr, int N);

// Print abundances present in network.
int print_abundances(const struct tnn* network);

// Print result data.
int print_results(const struct rate_library* rates, const struct tnn* network,
                  const struct problem_parameters* params);

// Clean up all options given passed to code, freeing pointers, etc.
int options_clean(struct option_values options);

// Clean all simulation_properties pointers, close files, etc.
int simulation_properties_clean(struct simulation_properties sim_prop);

// Return simulation_properties from toml_result_t.
// AND: reads config_toml values to option_values.
struct simulation_properties
simulation_properties_create(toml_result_t simulation_toml,
                             toml_result_t config_toml,
                             struct option_values* opts);

// Validate information in simulation_properties, returning
// EXIT_SUCCESS(FAILURE).
int simulation_properties_validate(struct simulation_properties* sim_prop);

#endif
