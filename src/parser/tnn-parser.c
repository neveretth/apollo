#include "tnn-parser.h"

#include "validate.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <H5Dpublic.h>
#include <hdf5.h>
#include <hdf5_hl.h>

#define PF_ALLOC_CHUNK 24

struct tnn* network_create(struct simulation_properties sim_prop) {

    if (validate_file(sim_prop.network_file) == EXIT_FAILURE) {
        return NULL;
    }

    int size = 0;
    fscanf(sim_prop.network_file, "%i\n", &size);

    struct tnn* network = malloc(sizeof(struct tnn));

    network->f = malloc(sizeof(struct tnn_f));
    network->info = calloc(1, sizeof(struct tnn_info));
    network->fptr = malloc(sizeof(struct tnn_fptr));
    network->iptr = malloc(sizeof(struct tnn_iptr));

    network->iptr->z = malloc(size * sizeof(int));
    network->iptr->n = malloc(size * sizeof(int));

    network->fptr->aa = malloc(size * sizeof(float));
    network->fptr->x = malloc(size * sizeof(float));
    network->fptr->y = malloc(size * sizeof(float));
    network->fptr->mass_excess = malloc(size * sizeof(float));

    network->info->part_func = malloc(size * sizeof(float*));
    network->info->iso_label = malloc(size * sizeof(char*));

    char line[60];
    char iso_symbol[5];
    int z, n, a;
    float y, mass;
    float pf0, pf1, pf2, pf3, pf4, pf5, pf6, pf7;

    // Genuinely atrocious... do NOT let this into release.
    float temp[24] = {0.1f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f,
                      0.8f, 0.9f,  1.0f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f,
                      4.0f, 4.5f,  5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f};
    network->info->part_func_temp = malloc(sizeof(float) * PF_ALLOC_CHUNK);
    for (int i = 0; i < 24; i++) {
        network->info->part_func_temp[i] = temp[i];
    }

    int iso_idx = -1;
    int iso_subindex = 3;

    // if (sim_prop.verbose) {
    //     printf("\nData read in:\n");
    // }

    // Assume lines can contain up to 60 characters.
    while (fgets(line, 60, sim_prop.network_file)) {
        iso_subindex++;
        if (iso_subindex == 4) {
            iso_subindex = 0;
            iso_idx++;
            // Scan and parse a title line
            sscanf(line, "%s %d %d %d %f %f", iso_symbol, &a, &z, &n, &y,
                   &mass);
            // if (sim_prop.verbose) {
            //     printf("\n%s %d %d %d %f %f\n", iso_symbol, a, z, n, y, mass);
            // }
            network->iptr->z[iso_idx] = z;
            network->iptr->n[iso_idx] = n;
            network->fptr->aa[iso_idx] = (float)a;
            network->fptr->y[iso_idx] = y;
            network->fptr->x[iso_idx] =
                network->fptr->aa[iso_idx] * network->fptr->y[iso_idx];
            network->fptr->mass_excess[iso_idx] = mass;
            network->info->iso_label[iso_idx] =
                malloc(sizeof(char) * (strlen(iso_symbol) + 1));
            strcpy(network->info->iso_label[iso_idx], iso_symbol);
        } else {
            // Scan and parse a partition function line.
            sscanf(line, "%f %f %f %f %f %f %f %f", &pf0, &pf1, &pf2, &pf3,
                   &pf4, &pf5, &pf6, &pf7);
            // if (sim_prop.verbose) {
            //     printf("%f %f %f %f %f %f %f %f\n", pf0, pf1, pf2, pf3, pf4,
            //            pf5, pf6, pf7);
            // }
            int tin = iso_subindex - 1;
            if (network->info->part_func[iso_idx] == NULL) {
                network->info->part_func[iso_idx] =
                    malloc(sizeof(float) * PF_ALLOC_CHUNK);
            }
            network->info->part_func[iso_idx][8 * (tin)] = pf0;
            network->info->part_func[iso_idx][8 * (tin) + 1] = pf1;
            network->info->part_func[iso_idx][8 * (tin) + 2] = pf2;
            network->info->part_func[iso_idx][8 * (tin) + 3] = pf3;
            network->info->part_func[iso_idx][8 * (tin) + 4] = pf4;
            network->info->part_func[iso_idx][8 * (tin) + 5] = pf5;
            network->info->part_func[iso_idx][8 * (tin) + 6] = pf6;
            network->info->part_func[iso_idx][8 * (tin) + 7] = pf7;
        }

        network->info->number_species = iso_idx + 1;
    }

    return network;
}
