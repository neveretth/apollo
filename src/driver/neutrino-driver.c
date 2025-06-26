#include "neutrino-driver.h"

#include "../parse-data.h"
#include "kernel-driver.h"

#include <stdlib.h>
#include <math.h>

int neutrino_debug(struct option_values options) {
    struct neunet* network = neunet_create(options);

    network->f->cycle_max = (int)1e9;
    network->f->t_end = 1.0e-00;
    network->f->t_w0 = 1.0e-10;
    network->f->EpsA = 1.0e-06;
    network->f->EpsR = 1.0e-04;
    network->f->t = 0.0;
    network->f->dt = 1.0e-15;
    network->f->dt_new = network->f->dt;
    network->f->tol = 0.0;
    network->f->g_a = 7.5e-01;
    network->f->g_b = 1.0e+02;
    network->f->g_c = sqrt(50.0);
    network->f->err = 0.0;
    options.halt = 10000;

    if (neunet_integrate_network(network, options) == EXIT_FAILURE) {
        neunet_destroy(&network);
        return EXIT_FAILURE;
    }

    neunet_print(network);

    neunet_destroy(&network);
    
    return EXIT_SUCCESS;
}
