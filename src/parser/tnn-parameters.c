#include "tnn-parameters.h"

#include <stdio.h>
#include <stdlib.h>

struct problem_parameters*
problem_parameters_create(struct rate_library* rates, struct tnn* network,
                          struct option_values options) {
    struct problem_parameters* params =
        calloc(1, sizeof(struct problem_parameters));

    // Find for each isotope all reactions that change its population.  This
    // analysis of the network is required only once at the very beginning of
    // the calculation (provided that the network species and reactions remain
    // the same for the entire calculation). The work is done by the function
    // reaction_mask_create().

    // Number of F+ and F- components for each isotope
    params->f_plus_number = malloc(sizeof(int) * network->info->number_species);
    params->f_minus_number =
        malloc(sizeof(int) * network->info->number_species);

    int* temp_int1 =
        calloc(network->info->number_species * rates->number_reactions / 2,
               sizeof(int));
    int* temp_int2 =
        calloc(network->info->number_species * rates->number_reactions / 2,
               sizeof(int));

    // reaction_mask_create() was formerly parse_f()
    params->reaction_mask = reaction_mask_create(rates, network, params,
                                                 options, temp_int1, temp_int2);

    // Create 1D arrays to hold non-zero F+ and F- for all reactions for all
    // isotopes, the arrays holding the species factors params->f_plus_factor
    // and params->f_minus_factor, and also arrays to hold their sums for each
    // isotope. Note that parseF() must be run first because it determines
    // params->f_plus_total and params->f_minus_total.

    params->f_plus = malloc(sizeof(float) * params->f_plus_total);
    params->f_minus = malloc(sizeof(float) * params->f_minus_total);
    params->f_plus_factor = malloc(sizeof(float) * params->f_plus_total);
    params->f_minus_factor = malloc(sizeof(float) * params->f_minus_total);
    params->f_plus_sum = malloc(sizeof(float) * network->info->number_species);
    params->f_minus_sum = malloc(sizeof(float) * network->info->number_species);

    // Arrays that hold the index of the boundary between different isotopes in
    // the f_plus and f_minus 1D arrays. Since f_plus_max and f_plus_min are
    // related, and likewise f_minus_max and f_minus_min are related, we will
    // only need to pass f_plus_max and f_minus_max to the kernel.

    params->f_plus_max = malloc(sizeof(int) * network->info->number_species);
    params->f_plus_min = malloc(sizeof(int) * network->info->number_species);
    params->f_minus_max = malloc(sizeof(int) * network->info->number_species);
    params->f_minus_min = malloc(sizeof(int) * network->info->number_species);

    // Create 1D arrays that will be used to map finite F+ and F- to the Flux
    // array.

    params->f_plus_isotope_cut =
        malloc(sizeof(int) * network->info->number_species);
    params->f_minus_isotope_cut =
        malloc(sizeof(int) * network->info->number_species);

    params->f_plus_isotope_idx =
        (int*)malloc(sizeof(int) * params->f_plus_total);
    params->f_minus_isotope_idx =
        (int*)malloc(sizeof(int) * params->f_minus_total);

    // Create 1D arrays that will hold the index of the isotope for the F+ or F-
    // term
    params->f_plus_map = malloc(sizeof(int) * params->f_plus_total);
    params->f_minus_map = malloc(sizeof(int) * params->f_minus_total);

    params->f_plus_isotope_cut[0] = params->f_plus_number[0];
    params->f_minus_isotope_cut[0] = params->f_minus_number[0];
    for (int i = 1; i < network->info->number_species; i++) {
        params->f_plus_isotope_cut[i] =
            params->f_plus_number[i] + params->f_plus_isotope_cut[i - 1];
        params->f_minus_isotope_cut[i] =
            params->f_minus_number[i] + params->f_minus_isotope_cut[i - 1];
    }

    int current_iso = 0;
    for (int i = 0; i < params->f_plus_total; i++) {
        params->f_plus_isotope_idx[i] = current_iso;
        if (i == (params->f_plus_isotope_cut[current_iso] - 1)) {
            current_iso++;
        }
    }

    current_iso = 0;
    for (int i = 0; i < params->f_minus_total; i++) {
        params->f_minus_isotope_idx[i] = current_iso;
        if (i == (params->f_minus_isotope_cut[current_iso] - 1))
            current_iso++;
    }

    // Diagnostic output
    if (options.verbose) {
        printf("\n\n\nMAX F+ and F- INDEX FOR EACH ISOTOPE:\n");
        for (int i = 0; i < network->info->number_species; i++) {
            printf("\nIsotope index = %d  %s  Max index F+ = %d  Max index F- "
                   "= %d",
                   i, network->info->iso_label[i],
                   params->f_plus_isotope_cut[i] - 1,
                   params->f_minus_isotope_cut[i] - 1);
        }
    }

    for (int i = 0; i < params->f_plus_total; i++) {
        params->f_plus_map[i] = temp_int1[i];
    }

    for (int i = 0; i < params->f_minus_total; i++) {
        params->f_minus_map[i] = temp_int2[i];
    }

    // Populate the params->f_plus_min and params->f_plus_max arrays
    params->f_plus_min[0] = 0;
    params->f_plus_max[0] = params->f_plus_number[0] - 1;
    for (int i = 1; i < network->info->number_species; i++) {
        params->f_plus_min[i] = params->f_plus_max[i - 1] + 1;
        params->f_plus_max[i] =
            params->f_plus_min[i] + params->f_plus_number[i] - 1;
    }
    // Populate the params->f_minus_min and params->f_minus_max arrays
    params->f_minus_min[0] = 0;
    params->f_minus_max[0] = params->f_minus_number[0] - 1;
    for (int i = 1; i < network->info->number_species; i++) {
        params->f_minus_min[i] = params->f_minus_max[i - 1] + 1;
        params->f_minus_max[i] =
            params->f_minus_min[i] + params->f_minus_number[i] - 1;
    }

    // Populate the params->f_plus_factor and params->f_minus_factor arrays that
    // hold the factors counting the number of occurences of the species in the
    // reaction.  Note that this can only be done after parseF() has been run to
    // give params->reaction_mask[i][j].

    int temp_count_plus = 0;
    int temp_count_minus = 0;
    for (int i = 0; i < network->info->number_species; i++) {
        for (int j = 0; j < rates->number_reactions; j++) {
            if (params->reaction_mask[i][j] > 0) {
                params->f_plus_factor[temp_count_plus] =
                    (float)params->reaction_mask[i][j];
                // 	printf("\n F+ temp_count_plus=%d
                // i=%d j=%d params->f_plus_factor=%3.1f", temp_count_plus,
                // i, j, params->f_plus_factor[temp_count_plus]);
                temp_count_plus++;
            } else if (params->reaction_mask[i][j] < 0) {
                params->f_minus_factor[temp_count_minus] =
                    -(float)params->reaction_mask[i][j];
                // 	printf("\n F-
                // temp_count_minus=%d i=%d j=%d params->f_minus_factor=%3.1f",
                // temp_count_minus, 					   i, j,
                // params->f_minus_factor[temp_count_minus]);
                temp_count_minus++;
            }
        }
    }

    if (options.verbose) {
        printf("\n\n\n---------- %d NON-VANISHING F+ SOURCE TERMS ----------\n",
               params->f_plus_total);
        printf("\ndY[%s]/dt = dY[%d]/dt F+ source terms (%d):",
               network->info->iso_label[params->f_plus_isotope_idx[0]],
               params->f_plus_isotope_idx[0],
               params->f_plus_number[params->f_plus_isotope_idx[0]]);
        for (int i = 0; i < params->f_plus_total; i++) {
            printf("\n   Isotope index = %d F+ index = %d Reac index = %d  %s",
                   params->f_plus_isotope_idx[i], i, params->f_plus_map[i],
                   rates->reaction_label[params->f_plus_map[i]]);
            if (i ==
                    (params->f_plus_isotope_cut[params->f_plus_isotope_idx[i]] -
                     1) &&
                i != params->f_plus_total - 1) {
                printf("\n");
                printf(
                    "\ndY[%s]/dt = dY[%d]/dt F+ source terms (%d):",
                    network->info->iso_label[params->f_plus_isotope_idx[i + 1]],
                    params->f_plus_isotope_idx[i + 1],
                    params->f_plus_number[params->f_plus_isotope_idx[i + 1]]);
            }
        }

        printf("\n\n\n---------- %d NON-VANISHING F- SOURCE TERMS ----------\n",
               params->f_minus_total);
        printf("\ndY[%s]/dt = dY[%d]/dt F- source terms (%d):",
               network->info->iso_label[params->f_minus_isotope_idx[0]],
               params->f_minus_isotope_idx[0],
               params->f_minus_number[params->f_minus_isotope_idx[0]]);
        for (int i = 0; i < params->f_minus_total; i++) {
            printf("\n   Isotope index = %d F- index = %d Reac index=%d  %s",
                   params->f_minus_isotope_idx[i], i, params->f_minus_map[i],
                   rates->reaction_label[params->f_minus_map[i]]);
            if (i ==
                    (params
                         ->f_minus_isotope_cut[params->f_minus_isotope_idx[i]] -
                     1) &&
                i != params->f_minus_total - 1) {
                printf("\n");
                printf(
                    "\ndY[%s]/dt = dY[%d]/dt F- source terms (%d):",
                    network->info
                        ->iso_label[params->f_minus_isotope_idx[i + 1]],
                    params->f_minus_isotope_idx[i + 1],
                    params->f_minus_number[params->f_minus_isotope_idx[i + 1]]);
            }
        }

        printf("\n");
    }

    free(temp_int1);
    free(temp_int2);

    return params;
}

int** reaction_mask_create(struct rate_library* rates, struct tnn* network,
                           struct problem_parameters* params,
                           struct option_values options, int* temp_int1,
                           int* temp_int2) {
    int** mask = malloc(sizeof(int*) * network->info->number_species);
    for (int i = 0; i < network->info->number_species; i++) {
        mask[i] = malloc(sizeof(int) * rates->number_reactions);
    }

    if (options.verbose) {
        printf("\nUse parseF() to find F+ and F- flux components for each "
               "species:");
    }

    int increment_plus = 0;
    int increment_minus = 0;

    for (int i = 0; i < network->info->number_species; i++) {
        int total = 0;
        int f_plus_number = 0;
        int f_minus_number = 0;
        if (options.verbose) {
            printf("\n");
        }

        // Loop over all possible reactions for this isotope, finding those that
        // change its population up (contributing to F+) or down (contributing
        // to F-).
        for (int j = 0; j < rates->number_reactions; j++) {
            int l_total = 0;
            int r_total = 0;

            // Loop over reactants for this reaction
            for (int k = 0; k < rates->num_react_species[j]; k++) {
                if (network->iptr->z[i] == rates->reactant_z[j][k] &&
                    network->iptr->n[i] == rates->reactant_n[j][k])
                    l_total++;
            }

            // Loop over products for this reaction
            for (int k = 0; k < rates->num_products[j]; k++) {
                if (network->iptr->z[i] == rates->product_z[j][k] &&
                    network->iptr->n[i] == rates->product_n[j][k])
                    r_total++;
            }

            total = l_total - r_total;

            if (total > 0) { // Contributes to F- for this isotope
                f_minus_number++;
                mask[i][j] = -total;
                temp_int2[increment_minus + f_minus_number - 1] = j;
                if (options.verbose) {
                    printf("\n%s reacIndex=%d %s nReac=%d nProd=%d totL=%d "
                           "totR=%d tot=%d F-",
                           network->info->iso_label[i], j,
                           rates->reaction_label[j],
                           rates->num_react_species[j], rates->num_products[j],
                           l_total, r_total, total);
                }
            } else if (total < 0) { // Contributes to F+ for this isotope
                f_plus_number++;
                mask[i][j] = -total;
                temp_int1[increment_plus + f_plus_number - 1] = j;
                if (options.verbose) {
                    printf("\n%s reacIndex=%d %s nReac=%d nProd=%d totL=%d "
                           "totR=%d tot=%d F+",
                           network->info->iso_label[i], j,
                           rates->reaction_label[j],
                           rates->num_react_species[j], rates->num_products[j],
                           l_total, r_total, total);
                }
            } else { // Does not contribute to flux for this isotope
                mask[i][j] = 0;
            }
        }

        // Keep track of the total number of F+ and F- terms in the network for
        // all isotopes
        params->f_plus_total += f_plus_number;
        params->f_minus_total += f_minus_number;

        params->f_plus_number[i] = f_plus_number;
        params->f_minus_number[i] = f_minus_number;

        increment_plus += f_plus_number;
        increment_minus += f_minus_number;

        if (options.verbose) {
            printf("\n%d %s numF+ = %d numF- = %d", i,
                   network->info->iso_label[i], f_plus_number, f_minus_number);
        }
    }

    // Display some cases
    if (options.verbose) {
        printf("\n\nPART OF FLUX-ISOTOPE COMPONENT ARRAY (-n --> F-; +n --> F+ "
               "for "
               "given isotope):");
    }

    if (network->info->number_species != 16 &&
        network->info->number_species > 25) {
        // Comment out this part of the if-block for alpha network to prevent
        // warnings about index being out of bounds. (Doesn't matter in
        // calculation since this block is not reached if it is an alpha network
        // with 16 species, but the compile generates a long string of
        // warnings.)  Uncomment to show up to 26 species for larger networks.

        // 		printf("\n\nIndex
        // Reaction%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s",
        // 			network->info->iso_label[0],
        // network->info->iso_label[1], network->info->iso_label[2],
        // network->info->iso_label[3], network->info->iso_label[4],
        // network->info->iso_label[5], network->info->iso_label[6],
        // network->info->iso_label[7], network->info->iso_label[8],
        // network->info->iso_label[9], network->info->iso_label[10],
        // network->info->iso_label[11], network->info->iso_label[12],
        // network->info->iso_label[13],
        // network->info->iso_label[14], network->info->iso_label[15],
        // network->info->iso_label[16], network->info->iso_label[17],
        // network->info->iso_label[18], network->info->iso_label[19],
        // network->info->iso_label[20], network->info->iso_label[21],
        // network->info->iso_label[22], network->info->iso_label[23],
        // network->info->iso_label[24], network->info->iso_label[25]
        // 		);
        // 		for(int j=0; j<rates->number_reactions; j++)
        // 		{
        //
        // 			printf( "\n %4d %22s %4d %4d %4d %4d %4d %4d %4d
        // %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d
        //  %4d %4d %4d %4d %4d %4d",
        // 				j, rates->reaction_label[j], mask[0][j],
        // mask[1][j], mask[2][j],
        // mask[3][j], mask[4][j], mask[5][j], mask[6][j],
        // mask[7][j], mask[8][j], mask[9][j], mask[10][j],
        // mask[11][j], mask[12][j], mask[13][j], mask[14][j],
        // mask[15][j], mask[16][j], mask[17][j], mask[18][j],
        // mask[19][j], mask[20][j], mask[21][j], mask[22][j],
        // mask[23][j], mask[24][j], mask[25][j]
        // 			);
        // 		}

    } else if (network->info->number_species > 15) // For alpha networks
    {
        if (options.verbose) {
            printf("\n\nIndex               "
                   "Reaction%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s",
                   network->info->iso_label[0], network->info->iso_label[1],
                   network->info->iso_label[2], network->info->iso_label[3],
                   network->info->iso_label[4], network->info->iso_label[5],
                   network->info->iso_label[6], network->info->iso_label[7],
                   network->info->iso_label[8], network->info->iso_label[9],
                   network->info->iso_label[10], network->info->iso_label[11],
                   network->info->iso_label[12], network->info->iso_label[13],
                   network->info->iso_label[14], network->info->iso_label[15]);
            for (int j = 0; j < rates->number_reactions; j++) {

                printf(
                    "\n %4d %22s %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d "
                    "%4d %4d %4d %4d %4d",
                    j, rates->reaction_label[j], mask[0][j], mask[1][j],
                    mask[2][j], mask[3][j], mask[4][j], mask[5][j], mask[6][j],
                    mask[7][j], mask[8][j], mask[9][j], mask[10][j],
                    mask[11][j], mask[12][j], mask[13][j], mask[14][j],
                    mask[15][j]);
            }
        }
    }

    if (options.verbose) {
        printf("\n\nFLUX SPARSENESS: Non-zero F+ = %d; Non-zero F- = %d, out "
               "of %d "
               "x %d = %d possibilities.\n",
               params->f_plus_total, params->f_minus_total,
               rates->number_reactions, network->info->number_species,
               rates->number_reactions * network->info->number_species);
    }

    return mask;
}

int problem_parameters_update(struct problem_parameters* params,
                              struct rate_library* rates, struct tnn* network) {
    
    params->f_plus_total = 0;
    params->f_minus_total = 0;

    int* temp_int1 =
        calloc(network->info->number_species * rates->number_reactions / 2,
               sizeof(int));
    int* temp_int2 =
        calloc(network->info->number_species * rates->number_reactions / 2,
               sizeof(int));

    reaction_mask_update(params->reaction_mask, rates, network, params,
                         temp_int1, temp_int2);

    params->f_plus_isotope_cut[0] = params->f_plus_number[0];
    params->f_minus_isotope_cut[0] = params->f_minus_number[0];
    for (int i = 1; i < network->info->number_species; i++) {
        params->f_plus_isotope_cut[i] =
            params->f_plus_number[i] + params->f_plus_isotope_cut[i - 1];
        params->f_minus_isotope_cut[i] =
            params->f_minus_number[i] + params->f_minus_isotope_cut[i - 1];
    }

    int current_iso = 0;
    for (int i = 0; i < params->f_plus_total; i++) {
        params->f_plus_isotope_idx[i] = current_iso;
        if (i == (params->f_plus_isotope_cut[current_iso] - 1)) {
            current_iso++;
        }
    }

    current_iso = 0;
    for (int i = 0; i < params->f_minus_total; i++) {
        params->f_minus_isotope_idx[i] = current_iso;
        if (i == (params->f_minus_isotope_cut[current_iso] - 1))
            current_iso++;
    }

    for (int i = 0; i < params->f_plus_total; i++) {
        params->f_plus_map[i] = temp_int1[i];
    }

    for (int i = 0; i < params->f_minus_total; i++) {
        params->f_minus_map[i] = temp_int2[i];
    }

    // Populate the params->f_plus_min and params->f_plus_max arrays
    params->f_plus_min[0] = 0;
    params->f_plus_max[0] = params->f_plus_number[0] - 1;
    for (int i = 1; i < network->info->number_species; i++) {
        params->f_plus_min[i] = params->f_plus_max[i - 1] + 1;
        params->f_plus_max[i] =
            params->f_plus_min[i] + params->f_plus_number[i] - 1;
    }
    // Populate the params->f_minus_min and params->f_minus_max arrays
    params->f_minus_min[0] = 0;
    params->f_minus_max[0] = params->f_minus_number[0] - 1;
    for (int i = 1; i < network->info->number_species; i++) {
        params->f_minus_min[i] = params->f_minus_max[i - 1] + 1;
        params->f_minus_max[i] =
            params->f_minus_min[i] + params->f_minus_number[i] - 1;
    }

    int temp_count_plus = 0;
    int temp_count_minus = 0;
    for (int i = 0; i < network->info->number_species; i++) {
        for (int j = 0; j < rates->number_reactions; j++) {
            if (params->reaction_mask[i][j] > 0) {
                params->f_plus_factor[temp_count_plus] =
                    (float)params->reaction_mask[i][j];
                temp_count_plus++;
            } else if (params->reaction_mask[i][j] < 0) {
                params->f_minus_factor[temp_count_minus] =
                    -(float)params->reaction_mask[i][j];
                temp_count_minus++;
            }
        }
    }

    free(temp_int1);
    free(temp_int2);

    return EXIT_SUCCESS;
}

void reaction_mask_update(int** mask, struct rate_library* rates,
                          struct tnn* network,
                          struct problem_parameters* params, int* temp_int1,
                          int* temp_int2) {

    int increment_plus = 0;
    int increment_minus = 0;

    for (int i = 0; i < network->info->number_species; i++) {
        int total = 0;
        int f_plus_number = 0;
        int f_minus_number = 0;

        for (int j = 0; j < rates->number_reactions; j++) {
            int l_total = 0;
            int r_total = 0;

            for (int k = 0; k < rates->num_react_species[j]; k++) {
                if (network->iptr->z[i] == rates->reactant_z[j][k] &&
                    network->iptr->n[i] == rates->reactant_n[j][k])
                    l_total++;
            }

            for (int k = 0; k < rates->num_products[j]; k++) {
                if (network->iptr->z[i] == rates->product_z[j][k] &&
                    network->iptr->n[i] == rates->product_n[j][k])
                    r_total++;
            }

            total = l_total - r_total;

            if (total > 0) {
                f_minus_number++;
                mask[i][j] = -total;
                temp_int2[increment_minus + f_minus_number - 1] = j;
            } else if (total < 0) { // Contributes to F+ for this isotope
                f_plus_number++;
                mask[i][j] = -total;
                temp_int1[increment_plus + f_plus_number - 1] = j;
            } else {
                mask[i][j] = 0;
            }
        }

        params->f_plus_total += f_plus_number;
        params->f_minus_total += f_minus_number;

        params->f_plus_number[i] = f_plus_number;
        params->f_minus_number[i] = f_minus_number;

        increment_plus += f_plus_number;
        increment_minus += f_minus_number;
    }
}
