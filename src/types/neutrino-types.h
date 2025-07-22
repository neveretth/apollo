#ifndef __NEUTRINO_TYPES_H
#define __NEUTRINO_TYPES_H

#include "../defs.h"

// Contain neunet real_t type data.
struct neunet_f {
    real_t mu;
    real_t kt;
    real_t rho;
    real_t ye;
    real_t cycle_max;
    real_t t_end;
    real_t t_w0;
    real_t EpsA;
    real_t EpsR;
    real_t t;
    real_t dt;
    real_t dt_new;
    real_t tol;
    real_t g_a;
    real_t g_b;
    real_t g_c;
    real_t err;
};

// Contain neunet real_t* type data.
struct neunet_fptr {
    real_t* ec;
    real_t* dv;
    real_t* n_eq;
    real_t* n_old;
};

// Contain neunet special data.
struct neunet_info {
    int num_groups;
    real_t** rate_in;
    real_t** rate_out;
};

// Contain data for neutrino network.
struct neunet {
    struct neunet_f* f;
    struct neunet_fptr* fptr;
    struct neunet_info* info;
};

// Free all data in neutrino network, setting source to NULL.
int neunet_destroy(struct neunet** __src);

// Print information in source.
int neunet_print(struct neunet* src);

// Return pointer to clone of source neunet.
struct neunet* neunet_clone(const struct neunet* src);

#endif
