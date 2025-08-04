#include "hydro-types.h"
#include "../types.h"

#include <stdlib.h>

int rt_hydro_mesh_destroy(struct rt_hydro_mesh** __src) {
    struct rt_hydro_mesh* src = *__src;
    for (int i = 0; i < src->dim[0]; i++) {
        for (int j = 0; j < src->dim[1]; j++) {
            free(src->temp[i][j]);
            free(src->density[i][j]);
        }
        free(src->temp[i]);
        free(src->density[i]);
    }
    free(src->temp);
    free(src->density);
    free(src);
    *__src = NULL;
    return EXIT_SUCCESS;
}
