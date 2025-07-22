#ifndef __NEUNET_PARSER_H
#define __NEUNET_PARSER_H

#include "../types.h"

// Return real_t** retrieved from HDF5 file (for neutrino ONLY)
real_t** hdf5_read_2d(const hid_t file, const char* dataset_name);

// Return real_t* retrieved from HDF5 file (for neutrino ONLY)
real_t* hdf5_read_1d(const hid_t file, const char* dataset_name);

// Return real_t retrieved from HDF5 file (for neutrino ONLY)
real_t hdf5_read(const hid_t file, const char* dataset_name);

// Return pointer to neutrino network (neunet) with the data in the file passed
// to the program.
struct neunet* neunet_create(struct simulation_properties sim_prop);

#endif
