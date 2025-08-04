#include "validate.h"

#include <stdlib.h>
#include <string.h>

int validate_file(FILE* file) {
    char* filesig = malloc(256 * sizeof(char));
    fscanf(file, "%s\n", filesig);
    if (strcmp(filesig, "%AAD") != 0) {
        printf("==apollo== ERROR: invalid header in file\n");
        return EXIT_FAILURE;
    }
    free(filesig);
    return EXIT_SUCCESS;
}
