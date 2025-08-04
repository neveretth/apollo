#ifndef __HYDRO_TYPES_H
#define __HYDRO_TYPES_H

#include "../defs.h"

// R3 (3D) hydro mesh
struct rt_hydro_mesh {
    real_t*** temp;
    real_t*** density;
    real_t volume;
    real_t dt;
    real_t t_end;
    real_t h;
    int dim[3];
};
#define hydro_mesh rt_hydro_mesh

int rt_hydro_mesh_destroy(struct rt_hydro_mesh** __src);

#endif
