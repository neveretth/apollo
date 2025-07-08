#ifndef __HYDRO_TYPES_H
#define __HYDRO_TYPES_H

struct linear_hydro_mesh {
    float* temp;
    float* density;
    float volume;
    float dt;
    float t_end;
    float h;
    int dim;
};
// This is to be deprecated in some later version...
#define hydro_mesh linear_hydro_mesh

struct flat_hydro_mesh {
    float** temp;
    float** density;
    float volume;
    float dt;
    float t_end;
    float h;
    int dim[2];
};

int hydro_mesh_destroy(struct hydro_mesh** __src);

void hydro_mesh_print(struct hydro_mesh* mesh);

void flat_hydro_mesh_print(struct flat_hydro_mesh* mesh);

int flat_hydro_mesh_destroy(struct flat_hydro_mesh** __src);

#endif
