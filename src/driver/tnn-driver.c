#include "tnn-driver.h"

#include "../parser.h"
#include "kernel-driver.h"

#include <stdlib.h>

int thermo_debug(struct option_values options) {
    struct rate_library* rates = rate_library_create(options);

    struct tnn* network = network_create(options);

    // ************************

    options.halt = 1000;

    network_print(network);

    // Set the temperature in units of 10^9 K and density in units of g/cm^3.
    // The temperature and density will be passed from the hydro code in an
    // operator-split coupling of this network to hydro. These will be used to
    // calculate the reaction rates in the network on the GPU. Since we are
    // assuming operator splitting, the temperature and density are assumed
    // constant for the entire network integration on the gPU.

    network->f->t9 = 6.0f;
    network->f->rho = 1.0e8;

    // Set the range of time integration and the initial timestep.  In an
    // operator-split coupling tmax will come from the hydro and dt_init will
    // likely be the last timestep of the previous network integration (for the
    // preceding hydro timestep).

    network->f->t_max = 1e-11;
    network->f->dt_init = 1e-17;

    float density[3];

    density[0] = 1.0f;
    density[1] = network->f->rho;
    density[2] = network->f->rho * network->f->rho;

    for (int i = 0; i < rates->number_reactions; i++) {
        rates->prefactor[i] *= (density[rates->num_react_species[i] - 1]);
    }

    // Find for each isotope all reactions that change its population.  This
    // analysis of the network is required only once at the very beginning of
    // the calculation (provided that the network species and reactions remain
    // the same for the entire calculation). The work is done by the function
    // parseF().

    struct problem_parameters* params =
        problem_parameters_create(rates, network, options);

    if (tnn_integrate_network(rates, network, params, options) == EXIT_FAILURE) {
        rate_library_destroy(&rates);
        network_destroy(&network);
        problem_parameters_destroy(&params, network->info->number_species);
        return EXIT_FAILURE;
    }

    rate_library_print(rates, network);

    print_results(rates, network, params);

    print_abundances(network);

    // ************************

    rate_library_destroy(&rates);
    problem_parameters_destroy(&params, network->info->number_species);
    network_destroy(&network);
    return EXIT_SUCCESS;
}
