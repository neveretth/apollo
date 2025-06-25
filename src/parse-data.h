#ifndef __PARSE_DATA_H
#define __PARSE_DATA_H

#include "types.h"

// Return pointer to rate_library struct containing rate library data. This is
// required only once, at the beginning of the entire calculation.
struct rate_library* rate_library_create(struct option_values options);

// Return pointer to network struct containing the component network structures.
// This is required only once at the very beginning of the calculation.
struct tnn* network_create(struct option_values options);

// Return float** retrieved from HDF5 file (for neutrino ONLY)
float** hdf5_read_2d(const hid_t file, const char* dataset_name);

// Return float* retrieved from HDF5 file (for neutrino ONLY)
float* hdf5_read_1d(const hid_t file, const char* dataset_name);

// Return float retrieved from HDF5 file (for neutrino ONLY)
float hdf5_read(const hid_t file, const char* dataset_name);

// Return pointer to neutrino network (neunet) with the data in the file passed
// to the program.
struct neunet* neunet_create(struct option_values options);

#endif
