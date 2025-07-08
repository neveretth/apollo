#include "args.h"
#include "driver.h"

#include <stdlib.h>

int main(int argc, char** argv) {
    // May exit here if certain conditions are not met...
    struct option_values options = parse_args(argc, argv);

    if (options.rocm_debug) {
        printf("==apollo== debug not implemented.");
        options_clean(options);
        return EXIT_SUCCESS;
    }

    if (options.thermo_debug) {
        if (thermo_debug(options) == EXIT_FAILURE) {
            options_clean(options);
            return EXIT_FAILURE;
        }
    }
    if (options.neutrino_debug) {
        if (neutrino_debug(options) == EXIT_FAILURE) {
            options_clean(options);
            return EXIT_FAILURE;
        }
    }
    if (options.hydro_debug) {
        if (hydro_debug(options) == EXIT_FAILURE) {
            options_clean(options);
            return EXIT_FAILURE;
        }
    }
    if (options.full) {
        if (full_simple(options) == EXIT_FAILURE) {
            options_clean(options);
            return EXIT_FAILURE;
        }
    }

    options_clean(options);

    return EXIT_SUCCESS;
}
