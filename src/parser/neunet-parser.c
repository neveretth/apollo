#include "neunet-parser.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

// "Temporary" workaround until HDF5 is fully integrated..
#define NUM_GROUPS 40

// Yet another "temporary" workaround.
#define scattering_kernel 1

float** hdf5_read_2d(const hid_t file, const char* dataset_name) {
    float* data = malloc(sizeof(float*) * NUM_GROUPS * NUM_GROUPS);
    hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dataset);
    float** datagrid = malloc(sizeof(float*) * NUM_GROUPS);
    for (int i = 0; i < NUM_GROUPS; i++) {
        datagrid[i] = malloc(sizeof(float) * NUM_GROUPS);
        for (int k = 0; k < NUM_GROUPS; k++) {
            datagrid[i][k] = data[(i * NUM_GROUPS) + k];
        }
    }
    return datagrid;
}

float* hdf5_read_1d(const hid_t file, const char* dataset_name) {
    float* data = malloc(sizeof(float) * NUM_GROUPS);
    hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dataset);
    return data;
}

float hdf5_read(const hid_t file, const char* dataset_name) {
    float data;
    hid_t dataset = H5Dopen2(file, dataset_name, H5P_DEFAULT);
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &data);
    H5Dclose(dataset);
    return data;
}

struct neunet* neunet_create(struct option_values options) {
    struct neunet* network = malloc(sizeof(struct neunet));
    network->info = malloc(sizeof(struct neunet_info));
    network->info->num_groups = NUM_GROUPS;
    network->f = malloc(sizeof(struct neunet_f));
    network->fptr = malloc(sizeof(struct neunet_fptr));

    network->fptr->ec =
        hdf5_read_1d(options.neutrino_file, "ProfileInfo/Energy");
    float* de = hdf5_read_1d(options.neutrino_file, "ProfileInfo/EnergyWidths");

    // Calculating dv based on ec and de
    network->fptr->dv = malloc(sizeof(float) * network->info->num_groups);
    for (int i = 0; i < network->info->num_groups; i++) {
        network->fptr->dv[i] = ((network->fptr->ec[i] + 0.5 * de[i]) *
                                    (network->fptr->ec[i] + 0.5 * de[i]) *
                                    (network->fptr->ec[i] + 0.5 * de[i]) -
                                (network->fptr->ec[i] - 0.5 * de[i]) *
                                    (network->fptr->ec[i] - 0.5 * de[i]) *
                                    (network->fptr->ec[i] - 0.5 * de[i])) /
                               3.0;
    }

    // Reading in Rho, T, Ye, Mu
    network->f->mu =
        hdf5_read(options.neutrino_file, "ProfileInfo/Chemical_Potential");
    network->f->kt =
        hdf5_read(options.neutrino_file, "ProfileInfo/Temperature");
    network->f->rho = hdf5_read(options.neutrino_file, "ProfileInfo/Density");
    network->f->ye =
        hdf5_read(options.neutrino_file, "ProfileInfo/Electron_Fraction");

    char dataset_name[128];

    if (scattering_kernel == 1) {
        strcpy(dataset_name,
               "Opacities/NES_Electron"); // Neutrino-Electron Scattering
    } else if (scattering_kernel == 2) {
        strcpy(dataset_name, "Opacities/NES_ElecAnti"); // AntiNeutrino-Electron
                                                        // Scattering
    } else if (scattering_kernel == 3) {
        strcpy(dataset_name,
               "Opacities/NES_MuTaAnti"); // MuonTau-AntiNeutrino Scattering
    } else if (scattering_kernel == 4) {
        strcpy(dataset_name, "Opacities/NES_MuonTau"); // Muon-Tau Scattering
    }

    network->info->rate_in = hdf5_read_2d(options.neutrino_file, dataset_name);

    // Invert the bastard. This is slow but it really doesn't matter (for now).
    // O(who cares)
    for (int i = 0; i < network->info->num_groups; i++) {
        for (int k = 1; k < network->info->num_groups / 2; k++) {
            float tmp = network->info->rate_in[i][k];
            network->info->rate_in[i][k] = network->info->rate_in[k][i];
            network->info->rate_in[k][i] = tmp;
        }
    }

    float c = 2.99792458e10;

    // Multiply each element by c
    for (size_t i = 0; i < network->info->num_groups; ++i) {
        for (size_t j = 0; j < network->info->num_groups; ++j) {
            network->info->rate_in[i][j] *= c;
        }
    }

    // Initialize n_eq based on mu and kt for the current model
    network->fptr->n_eq = malloc(sizeof(float) * network->info->num_groups);
    for (int i = 0; i < network->info->num_groups; i++) {
        network->fptr->n_eq[i] =
            1.0 /
            (exp((network->fptr->ec[i] - network->f->mu) / network->f->kt) +
             1.0);
    }

    for (int i = 0; i < network->info->num_groups; i++) {
        for (int j = 0; j < network->info->num_groups; j++) {
            if (j < i) {
                network->info->rate_in[i][j] =
                    network->info->rate_in[j][i] *
                    exp((network->fptr->ec[j] - network->fptr->ec[i]) /
                        network->f->kt);
            }
        }
    }

    network->info->rate_out =
        malloc(sizeof(float*) * network->info->num_groups);
    for (int i = 0; i < network->info->num_groups; i++) {
        network->info->rate_out[i] =
            malloc(sizeof(float) * network->info->num_groups);
    }

    for (int i = 0; i < network->info->num_groups; i++) {
        for (int j = 0; j < network->info->num_groups; j++) {
            network->info->rate_out[i][j] = network->info->rate_in[j][i];
        }
    }

    network->fptr->n_old = malloc(sizeof(float) * network->info->num_groups);

    return network;
}
