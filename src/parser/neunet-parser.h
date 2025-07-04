#ifndef __NEUNET_PARSER_H
#define __NEUNET_PARSER_H

#include "../types.h"

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
