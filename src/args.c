#include "args.h"

#include <H5Include.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* help_string =
    "apollo: An astrophysics solver.\n\n  help: apollo --help\n"
    "\n    Configuration of apollo is done in the simulation.toml file...\n"
    "    ...please refer to the Documentation for further information\n";


struct option_values parse_args(int argc, char** argv) {

    if (argc == 1) {
        printf("%s", help_string);
        exit(1);
    }

    struct option_values options;
    options.verbose = 0;
    options.rocm_accel = 0;
    options.config_file = NULL;
    options.simulation_file = NULL;

    // Not a lot of safety here...
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "--help") == 0) {
                printf("%s", help_string);
                exit(0);
            } else if (strcmp(argv[i], "-C") == 0) {
                options.config_file = argv[i + 1];
                i++;
            } else if (strcmp(argv[i], "-S") == 0) {
                options.simulation_file = argv[i + 1];
                i++;
            } else if (strcmp(argv[i], "--verbose") == 0) {
                options.verbose = 1;
            } else {
                printf("==apollo== ERROR: unknown flag: %s\n", argv[i]);
                printf("%s", help_string);
                exit(1);
            }
        }
    }

    int fail = 0;

    if (options.config_file == NULL) {
        printf("==apollo== error: no config file given \"-C\" <config.toml>\n");
        exit(1);
    }
    if (options.config_file == NULL) {
        printf("==apollo== error: no simulation file given \"-S\" "
               "<simulation.toml> \n");
        exit(1);
    }

    if (fail) {
        exit(1);
    }

    return options;
}
