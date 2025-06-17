#include "types.h"

#include <stdlib.h>

// Don't bother setting to NULL becuase they're all inaccessible anyways.
int rate_library_destroy(struct rate_library* src) {
    free(src->rg_class);
    free(src->rg_member_idx);
    free(src->reaction_library_class);
    free(src->num_react_species);
    free(src->num_products);
    free(src->is_ec);
    free(src->is_reverse);
    free(src->prefactor);
    free(src->q_value);
    free(src->rate);
    free(src->flux);
    free(src->p0);
    free(src->p1);
    free(src->p2);
    free(src->p3);
    free(src->p4);
    free(src->p5);
    free(src->p6);
    free(src->reactant_1);
    free(src->reactant_2);
    free(src->reactant_3);
    for (int i = 0; i < src->number_reactions; i++) {
        free(src->reactant_z[i]);
        free(src->reactant_n[i]);
        free(src->reactant_idx[i]);
        free(src->product_z[i]);
        free(src->product_n[i]);
        free(src->product_idx[i]);
        free(src->reaction_label[i]);
    }
    free(src->reactant_z);
    free(src->reactant_n);
    free(src->product_z);
    free(src->reactant_idx);
    free(src->product_n);
    free(src->product_idx);
    free(src->reaction_label);
    free(src);
    return EXIT_SUCCESS;
}

int network_destroy(struct network* src) {
    free(src->z);
    free(src->n);
    free(src->aa);
    free(src->y);
    free(src->x);
    free(src->mass_excess);
    for (int i = 0; i < src->number_species; i++) {
        free(src->iso_label[i]);
        free(src->part_func[i]);
    }
    free(src->iso_label);
    free(src->part_func);
    free(src->part_func_temp);
    free(src);
    return EXIT_SUCCESS;
}

int problem_parameters_destroy(struct problem_parameters* src, int reactions) {
    free(src->f_plus);
    free(src->f_minus);
    free(src->f_plus_factor);
    free(src->f_minus_factor);
    free(src->f_plus_sum);
    free(src->f_minus_sum);
    free(src->f_plus_max);
    free(src->f_plus_min);
    free(src->f_minus_max);
    free(src->f_minus_min);
    free(src->f_plus_isotope_cut);
    free(src->f_minus_isotope_cut);
    free(src->f_plus_number);
    free(src->f_minus_number);
    free(src->f_plus_isotope_idx);
    free(src->f_minus_isotope_idx);
    free(src->f_plus_map);
    free(src->f_minus_map);
    // f_plus_total and f_minus_total should match...
    for (int i = 0; i < reactions; i++) {
        free(src->reaction_mask[i]);
    }
    free(src->reaction_mask);
    free(src);
    return EXIT_SUCCESS;
}

int network_print(const struct network* network) {
    printf("\n\n%d ISOTOPES IN NETWORK:\n\n", network->number_species);
    printf(
        "Index  Isotope   A   Z   N  Abundance Y  MassFrac X  MassXS(MeV)\n");
    for (int i = 0; i < network->number_species; i++) {
        printf("%5d %8s %3d %3d %3d  %8.5e   %9.6f   %10.5f\n", i,
               network->iso_label[i], (int)network->aa[i], network->z[i],
               network->n[i], network->y[i], network->x[i],
               network->mass_excess[i]);
    }

    printf("\nPARTITION FUNCTION TABLE:\n");
    printf(
        "\n T9 = %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\
 %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f",
        network->part_func_temp[0], network->part_func_temp[1],
        network->part_func_temp[2], network->part_func_temp[3],
        network->part_func_temp[4], network->part_func_temp[5],
        network->part_func_temp[6], network->part_func_temp[7],
        network->part_func_temp[8], network->part_func_temp[9],
        network->part_func_temp[10], network->part_func_temp[11],
        network->part_func_temp[12], network->part_func_temp[13],
        network->part_func_temp[14], network->part_func_temp[15],
        network->part_func_temp[16], network->part_func_temp[17],
        network->part_func_temp[18], network->part_func_temp[19],
        network->part_func_temp[20], network->part_func_temp[21],
        network->part_func_temp[22], network->part_func_temp[23]);
    for (int j = 0; j < network->number_species; j++) {
        printf(
            "\n%-5s %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\
 %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f",
            network->iso_label[j], network->part_func[j][0],
            network->part_func[j][1], network->part_func[j][2],
            network->part_func[j][3], network->part_func[j][4],
            network->part_func[j][5], network->part_func[j][6],
            network->part_func[j][7], network->part_func[j][8],
            network->part_func[j][9], network->part_func[j][10],
            network->part_func[j][11], network->part_func[j][12],
            network->part_func[j][13], network->part_func[j][14],
            network->part_func[j][15], network->part_func[j][16],
            network->part_func[j][17], network->part_func[j][18],
            network->part_func[j][19], network->part_func[j][20],
            network->part_func[j][21], network->part_func[j][22],
            network->part_func[j][23]);
    }

    printf("\n");
    return EXIT_SUCCESS;
}

int rate_library_print(const struct rate_library* rates,
                       const struct network* network) {
    for (int i = 0; i < rates->number_reactions; i++) {
        printf(
            "%d %s Rate=%6.3e Y1=%6.3e Y2=%6.3e Y3=%6.3e Flux=%7.4e Q=%6.3f "
            "Prefac=%6.3e Reactants=%d\n",
            i, rates->reaction_label[i], rates->rate[i],
            network->y[rates->reactant_1[i]], network->y[rates->reactant_2[i]],
            network->y[rates->reactant_3[i]], rates->flux[i], rates->q_value[i],
            rates->prefactor[i], rates->num_react_species[i]);
    }
    return EXIT_SUCCESS;
}

int abundances_print(const struct network* network) {
    printf("\nIndex  Isotope   Abundance Y   Mass Frac X");
    float sum_x = 0.0f;
    for (int i = 0; i < network->number_species; i++) {
        float x = network->y[i] * network->aa[i];
        sum_x += x;
        printf("\n %4d     %4s    %8.4e    %8.4e", i, network->iso_label[i],
               network->y[i], x);
    }
    printf("\n\nsum X = %6.4f", sum_x);
    return EXIT_SUCCESS;
}
