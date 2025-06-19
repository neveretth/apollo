#include "parse-data.h"

#include <stdlib.h>
#include <string.h>

#define LABELSIZE 35
#define ALLOC_CHUNK 48 // for dynamic allocation, 48 is sufficient for alpha

struct rate_library* rate_library_create(struct option_values options) {

    struct rate_library* rates = malloc(sizeof(struct rate_library));

    rates->rg_member_idx = malloc(sizeof(int) * ALLOC_CHUNK);
    rates->rg_class = malloc(sizeof(int) * ALLOC_CHUNK);
    rates->reaction_label = malloc(sizeof(char*) * ALLOC_CHUNK);
    rates->reaction_library_class = malloc(sizeof(int) * ALLOC_CHUNK);
    rates->num_react_species = malloc(sizeof(int) * ALLOC_CHUNK);
    rates->num_products = malloc(sizeof(int) * ALLOC_CHUNK);
    rates->is_ec = malloc(sizeof(int) * ALLOC_CHUNK);
    rates->is_reverse = malloc(sizeof(int) * ALLOC_CHUNK);
    rates->prefactor = malloc(sizeof(float) * ALLOC_CHUNK);
    rates->q_value = malloc(sizeof(float) * ALLOC_CHUNK);
    rates->p0 = malloc(sizeof(float) * ALLOC_CHUNK);
    rates->p1 = malloc(sizeof(float) * ALLOC_CHUNK);
    rates->p2 = malloc(sizeof(float) * ALLOC_CHUNK);
    rates->p3 = malloc(sizeof(float) * ALLOC_CHUNK);
    rates->p4 = malloc(sizeof(float) * ALLOC_CHUNK);
    rates->p5 = malloc(sizeof(float) * ALLOC_CHUNK);
    rates->p6 = malloc(sizeof(float) * ALLOC_CHUNK);
    rates->reactant_n = malloc(sizeof(int*) * ALLOC_CHUNK);
    rates->reactant_z = malloc(sizeof(int*) * ALLOC_CHUNK);
    rates->product_n = malloc(sizeof(int*) * ALLOC_CHUNK);
    rates->product_z = malloc(sizeof(int*) * ALLOC_CHUNK);
    rates->reactant_1 = malloc(sizeof(int*) * ALLOC_CHUNK);
    rates->reactant_2 = malloc(sizeof(int*) * ALLOC_CHUNK);
    rates->reactant_3 = malloc(sizeof(int*) * ALLOC_CHUNK);
    rates->reactant_idx = malloc(sizeof(int*) * ALLOC_CHUNK);
    rates->product_idx = malloc(sizeof(int*) * ALLOC_CHUNK);

    rates->rate = malloc(sizeof(int) * ALLOC_CHUNK);
    rates->flux = malloc(sizeof(int) * ALLOC_CHUNK);

    char line[120];
    char reaction_token[LABELSIZE];
    float p0, p1, p2, p3, p4, p5, p6, q, sf;
    int i0, i1, i2, i3, i4, i5, i6;
    int ii[6];

    /*
    Read in the file line by line and parse into variables. The expected
    structure of each line is:
         float float float float float float float string
    ...each separated by a space, with no whitespace in the string.
    */

    int n = -1;
    int subindex = -1;

    if (options.verbose) {
        printf("\nData read in:\n\n");
    }

    while (fgets(line, 120, options.rate_library_file)) {
        subindex++;
        switch (subindex) {
        case 0:
            n++;
            sscanf(line, "%s %d %d %d %d %d %d %d %f %f", reaction_token, &i0,
                   &i1, &i2, &i3, &i4, &i5, &i6, &sf, &q);
            rates->reaction_label[n] =
                malloc(sizeof(char) * (strlen(reaction_token) + 1));
            strcpy(rates->reaction_label[n], reaction_token);
            rates->rg_class[n] = i0;
            rates->rg_member_idx[n] = i1;
            rates->reaction_library_class[n] = i2;
            rates->num_react_species[n] = i3;
            rates->num_products[n] = i4;
            rates->is_ec[n] = i5;
            rates->is_reverse[n] = i6;
            rates->prefactor[n] = sf;
            rates->q_value[n] = q;
            if (options.verbose) {
                printf("\n\nReaction Index = %d", n);
                printf("\nisReverseR = %d reaclibIndex = %d",
                       rates->is_reverse[n], rates->reaction_library_class[n]);
                printf("\n%s %d %d %d %d %d %d %d %f %f",
                       rates->reaction_label[n], rates->rg_class[n],
                       rates->rg_member_idx[n],
                       rates->reaction_library_class[n],
                       rates->num_react_species[n], rates->num_products[n],
                       rates->is_ec[n], rates->is_reverse[n],
                       rates->prefactor[n], rates->q_value[n]);
            }
            break;
        case 1:
            sscanf(line, "%f %f %f %f %f %f %f", &p0, &p1, &p2, &p3, &p4, &p5,
                   &p6);
            rates->p0[n] = p0;
            rates->p1[n] = p1;
            rates->p2[n] = p2;
            rates->p3[n] = p3;
            rates->p4[n] = p4;
            rates->p5[n] = p5;
            rates->p6[n] = p6;
            if (options.verbose) {
                printf("\n%f %f %f %f %f %f %f", rates->p0[n], rates->p1[n],
                       rates->p2[n], rates->p3[n], rates->p4[n], rates->p5[n],
                       rates->p6[n]);
            }
            break;
        case 2:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->reactant_z[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_react_species[n]; mm++) {
                rates->reactant_z[n][mm] = ii[mm];
                if (options.verbose) {
                    printf("\n  Reactant[%d]: Z=%d", mm,
                           rates->reactant_z[n][mm]);
                }
            }
            break;
        case 3:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->reactant_n[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_react_species[n]; mm++) {
                rates->reactant_n[n][mm] = ii[mm];
                if (options.verbose) {
                    printf("\n  Reactant[%d]: N=%d", mm,
                           rates->reactant_n[n][mm]);
                }
            }
            break;
        case 4:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->product_z[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_products[n]; mm++) {
                rates->product_z[n][mm] = ii[mm];
                if (options.verbose) {
                    printf("\n  Product[%d]: Z=%d", mm,
                           rates->product_z[n][mm]);
                }
            }
            break;
        case 5:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->product_n[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_products[n]; mm++) {
                rates->product_n[n][mm] = ii[mm];
                if (options.verbose) {
                    printf("\n  Product[%d]: N=%d", mm,
                           rates->product_n[n][mm]);
                }
            }
            break;
        case 6:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->reactant_idx[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_react_species[n]; mm++) {
                rates->reactant_idx[n][mm] = ii[mm];
                if (options.verbose) {
                    printf("\n  ReactantIndex[%d]: N=%d", mm,
                           rates->reactant_idx[n][mm]);
                }
            }
            break;
        case 7:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->product_idx[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_products[n]; mm++) {
                rates->product_idx[n][mm] = ii[mm];
                if (options.verbose) {
                    printf("\n  ProductIndex[%d]: N=%d", mm,
                           rates->product_idx[n][mm]);
                }
            }
            subindex = -1;
            break;
        }
    }
    rates->number_reactions = n + 1;

    for (int i = 0; i < rates->number_reactions; i++) {
        rates->reactant_1[i] = rates->reactant_idx[i][0];
        rates->reactant_2[i] = rates->reactant_idx[i][1];
        rates->reactant_3[i] = rates->reactant_idx[i][2];
    }
    if (options.verbose) {
        printf("\n");
    }

    return rates;
}

#define NETWORK_ALLOC_CHUNK 16
#define PF_ALLOC_CHUNK 24

struct thermo_network* network_create(struct option_values options) {

    struct thermo_network* network = malloc(sizeof(struct thermo_network));

    network->f = malloc(sizeof(struct tnn_f));
    network->info = malloc(sizeof(struct tnn_info));
    network->fptr = malloc(sizeof(struct tnn_fptr));
    network->iptr = malloc(sizeof(struct tnn_iptr));

    network->iptr->z = malloc(sizeof(int) * NETWORK_ALLOC_CHUNK);
    network->iptr->n = malloc(sizeof(int) * NETWORK_ALLOC_CHUNK);

    network->fptr->aa = malloc(sizeof(float) * NETWORK_ALLOC_CHUNK);
    network->fptr->x = malloc(sizeof(float) * NETWORK_ALLOC_CHUNK);
    network->fptr->y = malloc(sizeof(float) * NETWORK_ALLOC_CHUNK);
    network->fptr->mass_excess =
        malloc(sizeof(float) * NETWORK_ALLOC_CHUNK);

    network->info->part_func = malloc(sizeof(float*) * NETWORK_ALLOC_CHUNK);
    network->info->iso_label = malloc(sizeof(char*) * NETWORK_ALLOC_CHUNK);

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

    if (options.verbose) {
        printf("\nData read in:\n");
    }

    // Assume lines can contain up to 60 characters.
    while (fgets(line, 60, options.network_file)) {
        iso_subindex++;
        if (iso_subindex == 4) {
            iso_subindex = 0;
            iso_idx++;
            // Scan and parse a title line
            sscanf(line, "%s %d %d %d %f %f", iso_symbol, &a, &z, &n, &y,
                   &mass);
            if (options.verbose) {
                printf("\n%s %d %d %d %f %f\n", iso_symbol, a, z, n, y, mass);
            }
            network->iptr->z[iso_idx] = z;
            network->iptr->n[iso_idx] = n;
            network->fptr->aa[iso_idx] = (float)a;
            network->fptr->y[iso_idx] = y;
            network->fptr->x[iso_idx] = network->fptr->aa[iso_idx] *
                                              network->fptr->y[iso_idx];
            network->fptr->mass_excess[iso_idx] = mass;
            network->info->iso_label[iso_idx] =
                malloc(sizeof(char) * (strlen(iso_symbol) + 1));
            strcpy(network->info->iso_label[iso_idx], iso_symbol);
        } else {
            // Scan and parse a partition function line.
            sscanf(line, "%f %f %f %f %f %f %f %f", &pf0, &pf1, &pf2, &pf3,
                   &pf4, &pf5, &pf6, &pf7);
            if (options.verbose) {
                printf("%f %f %f %f %f %f %f %f\n", pf0, pf1, pf2, pf3, pf4,
                       pf5, pf6, pf7);
            }
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
