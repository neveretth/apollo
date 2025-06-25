#include "types.h"

#include <stdlib.h>
#include <H5Include.h>

int freenptr(void** ptr, int len) {
    if (ptr == NULL) {
        return EXIT_SUCCESS;
    }
    for (int i = 0; i < len; i++) {
        free(ptr[i]);
    }
    return EXIT_SUCCESS;
}

int print_abundances(const struct tnn* network) {
    printf("\n\nFINAL ABUNDANCES:\n");
    printf("\nIndex  Isotope   Abundance Y   Mass Frac X");
    float sum_x = 0.0f;
    for (int i = 0; i < network->info->number_species; i++) {
        float x = network->fptr->y[i] * network->fptr->aa[i];
        sum_x += x;
        printf("\n %4d     %4s    %8.4e    %8.4e", i,
               network->info->iso_label[i], network->fptr->y[i], x);
    }
    printf("\n\nsum X = %6.4f\n", sum_x);
    return EXIT_SUCCESS;
}

int print_results(const struct rate_library* rates, const struct tnn* network,
                  const struct problem_parameters* params) {
    printf("\n\nFINAL F+ VALUES:\n");
    for (int i = 0; i < params->f_plus_total; i++) {
        printf("\nF+[%d] = %7.4e  Increases Y[%s] through %s  MapIndex=%d  "
               "FplusFac=%3.1f",
               i, params->f_plus[i],
               network->info->iso_label[params->f_plus_isotope_idx[i]],
               rates->reaction_label[params->f_plus_map[i]],
               params->f_plus_map[i], params->f_plus_factor[i]);
    }
    printf("\n\n\nFINAL F- VALUES:\n");
    for (int i = 0; i < params->f_minus_total; i++) {
        printf("\nF-[%d] = %7.4e  Decreases Y[%s] through %s  MapIndex=%d  "
               "FminusFac=%3.1f",
               i, params->f_minus[i],
               network->info->iso_label[params->f_minus_isotope_idx[i]],
               rates->reaction_label[params->f_minus_map[i]],
               params->f_minus_map[i], params->f_minus_factor[i]);
    }

    printf("\n\n\nF+ and F- MIN AND MAX FOR EACH ISOTOPE:\n");
    for (int i = 0; i < network->info->number_species; i++) {
        printf("\n%3d %5s F+min=%3d F+max=%3d F-min=%3d F-max=%d", i,
               network->info->iso_label[i], params->f_plus_min[i],
               params->f_plus_max[i], params->f_minus_min[i],
               params->f_minus_max[i]);
    }

    printf("\n\n\nSUM OF FLUXES FOR EACH ISOTOPE:\n");

    float f_plus_total = 0.0f;
    float f_minus_total = 0.0f;
    for (int i = 0; i < network->info->number_species; i++) {
        printf("\n%3d %5s  sumF+=%10.4e  sumF-=%10.4e Fnet=%10.4e Y=%10.4e", i,
               network->info->iso_label[i], params->f_plus_sum[i],
               params->f_minus_sum[i],
               params->f_plus_sum[i] - params->f_minus_sum[i],
               network->fptr->y[i]);
        f_plus_total += params->f_plus_sum[i];
        f_minus_total += params->f_minus_sum[i];
    }

    printf("\n\ntotalF+ = %7.4e  totalF- = %7.4e", f_plus_total, f_minus_total);
    return EXIT_SUCCESS;
}

int options_clean(struct option_values options) {
    if (options.thermo_debug) {
        fclose(options.network_file);
        fclose(options.rate_library_file);
    }
    if (options.neutrino_debug) {
        H5Fclose(options.neutrino_file);
    }
    if (options.hydro_debug) {
    }
    return EXIT_SUCCESS;
}
