#ifndef __TYPES_H
#define __TYPES_H

#include "toml/tomlc17.h"
#include "types/hydro-types.h"
#include "types/neutrino-types.h"
#include "types/tnn-types.h"

#include <H5Include.h>
#include <stdio.h>

#include <stdbool.h>

// Contain values for all options when calling program.
struct option_values {
    FILE* rate_library_file;
    FILE* network_file;
    char* config_file;
    char* simulation_file;
    hid_t neutrino_file;
    int verbose;
    int rocm_accel;
    int rocm_debug;
    int halt;
    int thermo_debug;
    int neutrino_debug;
    int hydro_debug;
    int full;
};

struct simulation_properties {
    int sequence;
    int resolution[3];
    bool hydro;
    bool thermo;
    bool neutrino;
    FILE* rate_library_file;
    FILE* network_file;
    hid_t neutrino_file;
    bool output;
    float t_end;
    int output_tres;
    FILE* hydro_out_file;
    int (*hydro_temp_effect)(float*** data, int i_, int j_, int k_);
    int (*hydro_density_effect)(float*** data, int i_, int j_, int k_);
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

// Return simulation_properties from toml_result_t, exiting if failure occurs.
struct simulation_properties
simulation_properties_create(toml_result_t simulation_toml, toml_result_t config_toml);

// Validate information in simulation_properties, returning EXIT_SUCCESS(FAILURE).
int simulation_properties_validate(struct simulation_properties* sim_prop);

#endif
