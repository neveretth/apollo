#include "args.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* help_string = "apollo: An astrophysics solver.\n\n  usage: apollo "
                           "<rate-library-file> <network-file> --<flags>\n";

static char* help_extra_string = "\n  --help      print this info\n";

struct option_values parse_args(int argc, char** argv) {
    if (argc < 3) {
        if (argc > 1) {
            if (strcmp(argv[1], "--help") == 0) {
                printf("%s", help_string);
                printf("%s", help_extra_string);
                exit(1);
            }
        }
        printf("%s", help_string);
        exit(1);
    }

    struct option_values options;

    options.rate_library_file = fopen(argv[1], "r");
    if (options.rate_library_file == NULL) {
        printf("%s", help_string);
        exit(1);
    }
    printf("using rate library: %s\n", argv[1]);

    options.network_file = fopen(argv[2], "r");
    if (options.network_file == NULL) {
        printf("%s", help_string);
        exit(1);
    }
    printf("using network: %s\n", argv[2]);

    if (argc > 3) {
        for (int i = 3; i < argc; i++) {
            if (strcmp(argv[i], "--help") == 0) {
                printf("%s", help_string);
                printf("%s", help_extra_string);
                exit(1);
            } else {
                printf("Unknown flag: %s\n", argv[i]);
                printf("%s", help_string);
                exit(1);
            }
        }
    }

    return options;
}
