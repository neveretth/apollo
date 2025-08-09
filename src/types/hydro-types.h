#ifndef __HYDRO_TYPES_H
#define __HYDRO_TYPES_H

#include "../defs.h"

// R3 (3D) hydro mesh
// "block" ptrs are to the flat memory, the normal ptrs are to the ND data.
struct rt_hydro_mesh {
    real_t*** temp;
    real_t*** temp_block;
    real_t*** density;
    real_t*** density_block;
    real_t*** delta_temp;
    real_t*** delta_temp_block;
    real_t volume;
    real_t dt;
    real_t t_end;
    real_t h;
    int dim[3];
};
#define hydro_mesh rt_hydro_mesh

int rt_hydro_mesh_destroy(struct rt_hydro_mesh** __src);

#endif
