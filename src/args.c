#include "args.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <H5Include.h>

static char* help_string =
    "apollo: An astrophysics solver.\n\n  help: apollo --help\n";

static char* help_string_extra =
    "\n  --help            print this info\n"
    "  --rocm-accel      accelerate with rocm\n"
    "  --rocm-debug      debug rocm without running compute kernel\n"
    "  --full            run full-scale simulation\n";

static char* help_string_thermo =
    "\n *** Thermonuclear ***\n"
    "  --rate-library <file>   \n"
    "                    specify file to use for thermonuclear rate "
    "library\n"
    "  --network <file>  specify file to use for network intial state\n"
    "  --thermo-debug    debug thermonuclear computation\n";

static char* help_string_neutrino =
    "\n *** Neutrino ***\n"
    "  --neutrino-file <file>  \n"
    "                    specify file to use for initial neutrino data\n"
    "  --neutrino-debug  debug neutrino computation\n";

static char* help_string_hydro =
    "\n *** Hydro ***\n"
    "  --hydro-debug     debug hydro computation\n";

struct option_values parse_args(int argc, char** argv) {

    if (argc == 1) {
        printf("%s", help_string);
        exit(1);
    }

    struct option_values options;
    options.rate_library_file = NULL;
    options.network_file = NULL;
    options.neutrino_file = 0;
    options.verbose = 0;
    options.rocm_debug = 0;
    options.rocm_accel = 0;
    options.thermo_debug = 0;
    options.neutrino_debug = 0;
    options.hydro_debug = 0;
    options.full = 0;

    // Not a lot of safety here...
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "--help") == 0) {
                printf("%s", help_string);
                printf("%s", help_string_extra);
                printf("%s", help_string_thermo);
                printf("%s", help_string_neutrino);
                printf("%s", help_string_hydro);
                exit(0);
            } else if (strcmp(argv[i], "--verbose") == 0) {
                options.verbose = 1;
            } else if (strcmp(argv[i], "--rocm-debug") == 0) {
                options.rocm_debug = 1;
            } else if (strcmp(argv[i], "--rocm-accel") == 0) {
                options.rocm_accel = 1;
            } else if (strcmp(argv[i], "--thermo-debug") == 0) {
                options.thermo_debug = 1;
            } else if (strcmp(argv[i], "--neutrino-debug") == 0) {
                options.neutrino_debug = 1;
            } else if (strcmp(argv[i], "--hydro-debug") == 0) {
                options.hydro_debug = 1;
            } else if (strcmp(argv[i], "--full") == 0) {
                options.full = 1;
            } else if (strcmp(argv[i], "--rate-library") == 0) {
                options.rate_library_file = fopen(argv[i + 1], "r");
                if (options.rate_library_file == NULL) {
                    printf("%s", help_string);
                    exit(1);
                }
                printf("==apollo== using rate library: %s\n", argv[i + 1]);
                i++;
            } else if (strcmp(argv[i], "--network") == 0) {
                options.network_file = fopen(argv[i + 1], "r");
                if (options.network_file == NULL) {
                    printf("%s", help_string);
                    exit(1);
                }
                printf("==apollo== using network: %s\n", argv[i + 1]);
                i++;
            } else if (strcmp(argv[i], "--neutrino-file") == 0) {
                hid_t fapl;
                hsize_t size;
                if ((fapl = H5Pcreate(H5P_FILE_ACCESS)) == H5I_INVALID_HID) {
                    printf(
                        "==apollo== ERROR: cannot find file (dirty exit)");
                    exit(123);
                }
                if ((options.neutrino_file = H5Fopen(argv[i + 1], H5F_ACC_RDONLY, fapl)) ==
                    H5I_INVALID_HID) {
                    printf(
                        "==apollo== ERROR: cannot open file (dirty exit)");
                    exit(123);
                }
                printf("==apollo== using neutrino file: %s\n", argv[i + 1]);
                i++;
            } else {
                printf("==apollo== ERROR: unknown flag: %s\n", argv[i]);
                printf("%s", help_string);
                exit(1);
            }
        }
    }
    
    int fail = 0;

    // Review dependencies...
    if (options.thermo_debug) {
        if (options.rate_library_file == NULL || options.network_file == NULL) {
            printf("%s", help_string_thermo);
            fail++;
        }
    }
    if (options.neutrino_debug) {
        if (options.neutrino_file == 0) {
            printf("%s", help_string_neutrino);
            fail++;
        }
    }
    if (options.hydro_debug) {
    }
    if (options.full) {
        if (options.neutrino_file == 0) {
            printf("%s", help_string_neutrino);
            fail++;
        }
        if (options.rate_library_file == NULL || options.network_file == NULL) {
            printf("%s", help_string_thermo);
            fail++;
        }
    }

    if (fail) {
        exit(1);
    }

    return options;
}
