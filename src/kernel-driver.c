#include "kernel-driver.h"

#include "kernel.h"

#include <math.h>
#include <stdlib.h>

#define THIRD 0.3333333333333333

int integrate_network(struct rate_library* rates, struct network* network,
                      struct problem_parameters* params,
                      struct option_values options) {
    if (integration_kernel(
            rates->p0, rates->p1, rates->p2, rates->p3, rates->p4, rates->p5,
            rates->p6, rates->prefactor, rates->q_value, rates->rate,
            rates->flux, params->f_plus, params->f_minus, params->f_plus_factor,
            params->f_minus_factor, params->f_plus_sum, params->f_minus_sum,
            params->f_plus_max, params->f_minus_max, params->f_plus_map,
            params->f_minus_map, network->y, rates->num_react_species,
            rates->reactant_1, rates->reactant_2, rates->reactant_3,
            network->number_species, rates->number_reactions,
            params->f_plus_total, params->f_minus_total, options.t9,
            options.t_max, options.dt_init, options.halt) == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
