#include "rate-library-parser.h"

#include "validate.h"

#include <stdlib.h>
#include <string.h>

#define LABELSIZE 35

struct rate_library*
rate_library_create(struct simulation_properties sim_prop) {

    if (validate_file(sim_prop.rate_library_file) == EXIT_FAILURE) {
        return NULL;
    }

    // This _should_ probably be a library function.
    int size = 0;
    fscanf(sim_prop.rate_library_file, "%i\n", &size);

    struct rate_library* rates = malloc(sizeof(struct rate_library));

    rates->rg_member_idx = malloc(size * sizeof(int));
    rates->rg_class = malloc(size * sizeof(int));
    rates->reaction_label = malloc(size * sizeof(char*));
    rates->reaction_library_class = malloc(size * sizeof(int));
    rates->num_react_species = malloc(size * sizeof(int));
    rates->num_products = malloc(size * sizeof(int));
    rates->is_ec = malloc(size * sizeof(int));
    rates->is_reverse = malloc(size * sizeof(int));
    rates->prefactor = malloc(size * sizeof(real_t));
    rates->q_value = malloc(size * sizeof(real_t));
    rates->p0 = malloc(size * sizeof(real_t));
    rates->p1 = malloc(size * sizeof(real_t));
    rates->p2 = malloc(size * sizeof(real_t));
    rates->p3 = malloc(size * sizeof(real_t));
    rates->p4 = malloc(size * sizeof(real_t));
    rates->p5 = malloc(size * sizeof(real_t));
    rates->p6 = malloc(size * sizeof(real_t));
    rates->reactant_n = malloc(size * sizeof(int*));
    rates->reactant_z = malloc(size * sizeof(int*));
    rates->product_n = malloc(size * sizeof(int*));
    rates->product_z = malloc(size * sizeof(int*));
    rates->reactant_1 = malloc(size * sizeof(int*));
    rates->reactant_2 = malloc(size * sizeof(int*));
    rates->reactant_3 = malloc(size * sizeof(int*));
    rates->reactant_idx = malloc(size * sizeof(int*));
    rates->product_idx = malloc(size * sizeof(int*));

    rates->rate = malloc(size * sizeof(int));
    rates->flux = malloc(size * sizeof(int));

    char line[120];
    char reaction_token[LABELSIZE];
    real_t p0, p1, p2, p3, p4, p5, p6, q, sf;
    int i0, i1, i2, i3, i4, i5, i6;
    int ii[6];

    /*
    Read in the file line by line and parse into variables. The expected
    structure of each line is:
         real_t real_t real_t real_t real_t real_t real_t string
    ...each separated by a space, with no whitespace in the string.
    */

    int n = -1;
    int subindex = -1;

    // if (options.verbose) {
    //     printf("\nData read in:\n\n");
    // }

    while (fgets(line, 120, sim_prop.rate_library_file)) {
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
            // if (options.verbose) {
            //     printf("\n\nReaction Index = %d", n);
            //     printf("\nisReverseR = %d reaclibIndex = %d",
            //            rates->is_reverse[n],
            //            rates->reaction_library_class[n]);
            //     printf("\n%s %d %d %d %d %d %d %d %f %f",
            //            rates->reaction_label[n], rates->rg_class[n],
            //            rates->rg_member_idx[n],
            //            rates->reaction_library_class[n],
            //            rates->num_react_species[n], rates->num_products[n],
            //            rates->is_ec[n], rates->is_reverse[n],
            //            rates->prefactor[n], rates->q_value[n]);
            // }
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
            // if (options.verbose) {
            //     printf("\n%f %f %f %f %f %f %f", rates->p0[n], rates->p1[n],
            //            rates->p2[n], rates->p3[n], rates->p4[n],
            //            rates->p5[n], rates->p6[n]);
            // }
            break;
        case 2:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->reactant_z[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_react_species[n]; mm++) {
                rates->reactant_z[n][mm] = ii[mm];
                // if (options.verbose) {
                //     printf("\n  Reactant[%d]: Z=%d", mm,
                //            rates->reactant_z[n][mm]);
                // }
            }
            break;
        case 3:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->reactant_n[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_react_species[n]; mm++) {
                rates->reactant_n[n][mm] = ii[mm];
                // if (options.verbose) {
                //     printf("\n  Reactant[%d]: N=%d", mm,
                //            rates->reactant_n[n][mm]);
                // }
            }
            break;
        case 4:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->product_z[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_products[n]; mm++) {
                rates->product_z[n][mm] = ii[mm];
                // if (options.verbose) {
                //     printf("\n  Product[%d]: Z=%d", mm,
                //            rates->product_z[n][mm]);
                // }
            }
            break;
        case 5:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->product_n[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_products[n]; mm++) {
                rates->product_n[n][mm] = ii[mm];
                // if (options.verbose) {
                //     printf("\n  Product[%d]: N=%d", mm,
                //            rates->product_n[n][mm]);
                // }
            }
            break;
        case 6:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->reactant_idx[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_react_species[n]; mm++) {
                if (ii[mm] > 5) {
                    ii[mm] = 0;
                }
                rates->reactant_idx[n][mm] = ii[mm];
                // if (options.verbose) {
                    // printf("\n  ReactantIndex[%d]: N=%d", mm,
                    //        rates->reactant_idx[n][mm]);
                // }
            }
            break;
        case 7:
            sscanf(line, "%d %d %d %d", &ii[0], &ii[1], &ii[2], &ii[3]);
            rates->product_idx[n] = malloc(sizeof(int) * 4);
            for (int mm = 0; mm < rates->num_products[n]; mm++) {
                rates->product_idx[n][mm] = ii[mm];
                // if (options.verbose) {
                    // printf("\n  ProductIndex[%d]: N=%d", mm,
                    //        rates->product_idx[n][mm]);
                // }
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
    // if (options.verbose) {
    //     printf("\n");
    // }

    return rates;
}
