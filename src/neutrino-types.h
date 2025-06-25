#ifndef __NEUTRINO_TYPES_H
#define __NEUTRINO_TYPES_H

// Contain neunet float type data.
struct neunet_f {
    float mu;
    float kt;
    float rho;
    float ye;
    float cycle_max;
    float t_end;
    float t_w0;
    float EpsA;
    float EpsR;
    float t;
    float dt;
    float dt_new;
    float tol;
    float g_a;
    float g_b;
    float g_c;
    float err;
};

// Contain neunet float* type data.
struct neunet_fptr {
    float* ec;
    float* dv;
    float* n_eq;
    float* n_old;
};

// Contain neunet special data.
struct neunet_info {
    int num_groups;
    float** rate_in;
    float** rate_out;
};

// Contain data for neutrino network.
struct neunet {
    struct neunet_f* f;
    struct neunet_fptr* fptr;
    struct neunet_info* info;
};

int neunet_destroy(struct neunet** __src);

#endif
