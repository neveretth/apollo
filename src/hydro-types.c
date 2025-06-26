#include "hydro-types.h"

#include <stdio.h>
#include <stdlib.h>

int hydro_mesh_destroy(struct hydro_mesh** __src) {
    struct hydro_mesh* src = *__src;
    free(src->temp);
    free(src->density);
    free(src);
    *__src = NULL;
    return EXIT_SUCCESS;
}

void hydro_mesh_print(struct hydro_mesh* mesh) {
    for (int i = 0; i < mesh->dim; i++) {
        printf("%f ", mesh->temp[i]);
    }
    printf("\n");
}
