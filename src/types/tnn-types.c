#include "tnn-types.h"

#include "../types.h"

#include <stdlib.h>

int rate_library_destroy(struct rate_library** __src) {
    struct rate_library* src = *__src;
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
    freenptr((void**)src->reactant_z, src->number_reactions);
    free(src->reactant_z);
    freenptr((void**)src->reactant_n, src->number_reactions);
    free(src->reactant_n);
    freenptr((void**)src->reactant_idx, src->number_reactions);
    free(src->reactant_idx);
    freenptr((void**)src->product_z, src->number_reactions);
    free(src->product_z);
    freenptr((void**)src->product_idx, src->number_reactions);
    free(src->product_idx);
    freenptr((void**)src->product_n, src->number_reactions);
    free(src->product_n);
    freenptr((void**)src->reaction_label, src->number_reactions);
    free(src->reaction_label);
    free(src);
    *__src = NULL;
    return EXIT_SUCCESS;
}

int network_destroy(struct tnn** __src) {
    struct tnn* src = *__src;
    freenptr((void*)src->iptr, 2);
    freenptr((void*)src->fptr, 4);
    free(src->iptr);
    free(src->fptr);
    freenptr((void**)src->info->iso_label, src->info->number_species);
    free(src->info->iso_label);
    freenptr((void**)src->info->part_func, src->info->number_species);
    free(src->info->part_func);
    free(src->info->part_func_temp);
    free(src->info);
    free(src->f);
    free(src);
    *__src = NULL;
    return EXIT_SUCCESS;
}

int problem_parameters_destroy(struct problem_parameters** __src,
                               int reactions) {
    struct problem_parameters* src = *__src;
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
    freenptr((void**)src->reaction_mask, reactions);
    free(src->reaction_mask);
    free(src);
    *__src = NULL;
    return EXIT_SUCCESS;
}

int network_print(const struct tnn* network) {
    printf("\n\n%d ISOTOPES IN NETWORK:\n\n", network->info->number_species);
    printf(
        "Index  Isotope   A   Z   N  Abundance Y  MassFrac X  MassXS(MeV)\n");
    for (int i = 0; i < network->info->number_species; i++) {
        printf("%5d %8s %3d %3d %3d  %8.5e   %9.6f   %10.5f\n", i,
               network->info->iso_label[i], (int)network->fptr->aa[i],
               network->iptr->z[i], network->iptr->n[i], network->fptr->y[i],
               network->fptr->x[i], network->fptr->mass_excess[i]);
    }

    printf("\nPARTITION FUNCTION TABLE:\n");
    printf(
        "\n T9 = %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\
 %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f",
        network->info->part_func_temp[0], network->info->part_func_temp[1],
        network->info->part_func_temp[2], network->info->part_func_temp[3],
        network->info->part_func_temp[4], network->info->part_func_temp[5],
        network->info->part_func_temp[6], network->info->part_func_temp[7],
        network->info->part_func_temp[8], network->info->part_func_temp[9],
        network->info->part_func_temp[10], network->info->part_func_temp[11],
        network->info->part_func_temp[12], network->info->part_func_temp[13],
        network->info->part_func_temp[14], network->info->part_func_temp[15],
        network->info->part_func_temp[16], network->info->part_func_temp[17],
        network->info->part_func_temp[18], network->info->part_func_temp[19],
        network->info->part_func_temp[20], network->info->part_func_temp[21],
        network->info->part_func_temp[22], network->info->part_func_temp[23]);
    for (int j = 0; j < network->info->number_species; j++) {
        printf(
            "\n%-5s %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f\
 %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f %4.2f",
            network->info->iso_label[j], network->info->part_func[j][0],
            network->info->part_func[j][1], network->info->part_func[j][2],
            network->info->part_func[j][3], network->info->part_func[j][4],
            network->info->part_func[j][5], network->info->part_func[j][6],
            network->info->part_func[j][7], network->info->part_func[j][8],
            network->info->part_func[j][9], network->info->part_func[j][10],
            network->info->part_func[j][11], network->info->part_func[j][12],
            network->info->part_func[j][13], network->info->part_func[j][14],
            network->info->part_func[j][15], network->info->part_func[j][16],
            network->info->part_func[j][17], network->info->part_func[j][18],
            network->info->part_func[j][19], network->info->part_func[j][20],
            network->info->part_func[j][21], network->info->part_func[j][22],
            network->info->part_func[j][23]);
    }

    printf("\n");
    return EXIT_SUCCESS;
}

int rate_library_print(const struct rate_library* rates,
                       const struct tnn* network) {
    for (int i = 0; i < rates->number_reactions; i++) {
        printf("%d %s Rate=%6.3e Y1=%6.3e Y2=%6.3e Y3=%6.3e Flux=%7.4e Q=%6.3f "
               "Prefac=%6.3e Reactants=%d\n",
               i, rates->reaction_label[i], rates->rate[i],
               network->fptr->y[rates->reactant_1[i]],
               network->fptr->y[rates->reactant_2[i]],
               network->fptr->y[rates->reactant_3[i]], rates->flux[i],
               rates->q_value[i], rates->prefactor[i],
               rates->num_react_species[i]);
    }
    return EXIT_SUCCESS;
}

