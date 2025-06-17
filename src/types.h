#ifndef __TYPES_H
#define __TYPES_H

#include <stdio.h>

// Contain values for all options when calling program.
struct option_values {
    FILE* rate_library_file;
    FILE* network_file;
    int verbose;
};

// Contain rate library data.
// rg <-- means reaction group
struct rate_library {
    char** reaction_label;
    int* rg_class;
    int* rg_member_idx;
    int* reaction_library_class;
    int* num_react_species;
    int* num_products;
    int* is_ec; // electron capture
    int* is_reverse;
    float* prefactor; // statistical prefactor
    float* q_value;
    float* p0;
    float* p1;
    float* p2;
    float* p3;
    float* p4;
    float* p5;
    float* p6;
    int* reactant_1;
    int* reactant_2;
    int* reactant_3;
    int** reactant_z;
    int** reactant_n;
    int** product_z;
    int** product_n;
    int** reactant_idx;
    int** product_idx;
    int number_reactions;
};

struct network {
    int* z;
    int* n;
    float* aa;
    float* x;
    float* y;
    float* mass_excess;
    float** part_func;
    char** iso_label;
    int number_species;
};

int rate_library_destroy(struct rate_library* src);

int network_destroy(struct network* src);

#endif
