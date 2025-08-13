#include "hydro-types.h"
#include "../types.h"

#include <stdlib.h>

int rt_hydro_mesh_destroy(struct rt_hydro_mesh** __src) {
    struct rt_hydro_mesh* src = *__src;
    for (int k = 0; k < src->dim[ZDIM]; k++) {
        free(src->temp[k]);
        free(src->density[k]);
        free(src->delta_temp[k]);
        free(src->pressure[k]);
        free(src->velocity[k]);
        free(src->mean_mol_mass[k]);
    }
    free(src->temp);
    free(src->density);
    free(src->delta_temp);
    free(src->velocity);
    free(src->pressure);
    free(src->mean_mol_mass);

    free(src->temp_block);
    free(src->density_block);
    free(src->delta_temp_block);
    free(src->velocity_block);
    free(src->pressure_block);
    free(src->mean_mol_mass_block);

    free(src);
    *__src = NULL;
    return EXIT_SUCCESS;
}
