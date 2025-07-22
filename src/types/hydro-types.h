#ifndef __HYDRO_TYPES_H
#define __HYDRO_TYPES_H

#include "../defs.h"

struct linear_hydro_mesh {
    real_t* temp;
    real_t* density;
    real_t volume;
    real_t dt;
    real_t t_end;
    real_t h;
    int dim;
};
// This is to be deprecated in some later version...
#define hydro_mesh linear_hydro_mesh

struct flat_hydro_mesh {
    real_t** temp;
    real_t** density;
    real_t volume;
    real_t dt;
    real_t t_end;
    real_t h;
    int dim[2];
};

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

int hydro_mesh_destroy(struct hydro_mesh** __src);

void hydro_mesh_print(struct hydro_mesh* mesh);

void flat_hydro_mesh_print(struct flat_hydro_mesh* mesh);

int flat_hydro_mesh_destroy(struct flat_hydro_mesh** __src);

int rt_hydro_mesh_destroy(struct rt_hydro_mesh** __src);

#endif
