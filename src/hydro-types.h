#ifndef __HYDRO_TYPES_H
#define __HYDRO_TYPES_H

struct hydro_mesh {
    float* temp;
    float* density;
    float volume;
    float dt;
    float h;
    int dim;
};

int hydro_mesh_destroy(struct hydro_mesh** __src);

void hydro_mesh_print(struct hydro_mesh* mesh);

#endif
