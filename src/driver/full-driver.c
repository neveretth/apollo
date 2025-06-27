#include "full-driver.h"

#include "../parameters.h"
#include "../parse-data.h"
#include "kernel-driver.h"

#include <math.h>
#include <stdlib.h>

#define DIM 4

int full_simple(struct option_values options) {

    // Assume a linear arrangement of zones (1D) of size DIM

    struct tnn** thermo = malloc(DIM * sizeof(struct tnn*));
    struct rate_library** rates = malloc(DIM * sizeof(struct rate_library*));
    struct neunet** neutrino = malloc(DIM * sizeof(struct neunet*));
    for (int i = 0; i < DIM; i++) {
        thermo[i] = network_create(options);
        rates[i] = rate_library_create(options);
        neutrino[i] = neunet_create(options);
        // BLOODY AWFUL hack in lieu of a proper clone func...
        fseek(options.rate_library_file, 0, SEEK_SET);
        fseek(options.network_file, 0, SEEK_SET);
    }

    struct hydro_mesh* mesh = malloc(sizeof(struct hydro_mesh));
    mesh->dim = DIM;
    mesh->temp = malloc(mesh->dim * sizeof(float));
    mesh->density = malloc(mesh->dim * sizeof(float));

    // ************************

    mesh->dt = 1e-11;
    mesh->h = 10e+12;
    mesh->volume = 1;

    for (int i = 0; i < mesh->dim; i++) {
        neutrino[i]->f->cycle_max = (int)1e9;
        neutrino[i]->f->EpsA = 1.0e-06;
        neutrino[i]->f->EpsR = 1.0e-04;
        neutrino[i]->f->tol = 0.0;
        neutrino[i]->f->g_a = 7.5e-01;
        neutrino[i]->f->g_b = 1.0e+02;
        neutrino[i]->f->g_c = sqrt(50.0);
        neutrino[i]->f->dt = mesh->dt;
        neutrino[i]->f->dt_new = neutrino[i]->f->dt;
        thermo[i]->f->dt_init = 1e-17;
        // Density must be very low for the hydro to work correctly in 
        // it's current form... this is a bug.
        // thermo[i]->f->rho = 1.0e8;
        thermo[i]->f->rho = 1.0e0;
        // Apply some difference for the sake of testing.
        thermo[i]->f->t9 = 6.0f + (rand() % (2));
    }

    options.halt = 10000;

    float t = 0;
    float t_end = 1e-08;
    
    for (int i = 0; i < mesh->dim; i++) {
        print_abundances(thermo[i]);
        neunet_print(neutrino[i]);
    }

    // If this fails, it exits in a bloody massacre of memory leaks...
    struct problem_parameters* params;
    while (t < t_end) {
        printf("Time: %f/%f\n", t, t_end);
        // Transfer data from thermo/neutrino to hydro
        for (int i = 0; i < mesh->dim; i++) {
            // dt is not updated by the hydro timestep...
            //  but this is how the dataflow will work.
            thermo[i]->f->t_max = mesh->dt;
            neutrino[i]->f->t_end = mesh->dt;
            // The Thermonuclear code assumes a unit of 10^9, so account for
            // that. This is obviously inefficient as all hell, we will must fix
            // it before release!
            mesh->temp[i] = thermo[i]->f->t9 * 1e09;
            mesh->density[i] = thermo[i]->f->rho;
        }
        hydro_mesh_print(mesh);
        // Hydro integration step
        if (hydro_integrate_mesh(mesh, options) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
        hydro_mesh_print(mesh);
        // Transfer data from hydro to thermo/neutrino
        for (int i = 0; i < mesh->dim; i++) {
            thermo[i]->f->t9 = mesh->temp[i] / 1e09;
            thermo[i]->f->rho = mesh->density[i];
        }
        // Thermo/Neutrino integration step.
        for (int i = 0; i < mesh->dim; i++) {
            float density[3];
            density[0] = 1.0f;
            density[1] = thermo[i]->f->rho;
            density[2] = thermo[i]->f->rho * thermo[i]->f->rho;
            for (int k = 0; k < rates[i]->number_reactions; k++) {
                rates[i]->prefactor[k] *=
                    (density[rates[i]->num_react_species[k] - 1]);
            }
            params = problem_parameters_create(rates[i], thermo[i], options);
            if (tnn_integrate_network(rates[i], thermo[i], params, options) ==
                EXIT_FAILURE) {
                return EXIT_FAILURE;
            }
            if (neunet_integrate_network(neutrino[i], options) ==
                EXIT_FAILURE) {
                return EXIT_FAILURE;
            }
            problem_parameters_destroy(&params,
                                       thermo[i]->info->number_species);
        }
        t += mesh->dt;
    }
    printf("Finished!\n");

    for (int i = 0; i < mesh->dim; i++) {
        print_abundances(thermo[i]);
        neunet_print(neutrino[i]);
    }
    hydro_mesh_print(mesh);

    // ************************

    for (int i = 0; i < DIM; i++) {
        network_destroy(&(thermo[i]));
        rate_library_destroy(&(rates[i]));
        neunet_destroy(&(neutrino[i]));
    }
    hydro_mesh_destroy(&mesh);

    return EXIT_SUCCESS;
}
