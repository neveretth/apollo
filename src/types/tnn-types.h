#ifndef __TNN_TYPES_H
#define __TNN_TYPES_H

#include "../defs.h"

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
    real_t* rate;
    real_t* flux;
    real_t* prefactor; // statistical prefactor
    real_t* q_value;
    real_t* p0;
    real_t* p1;
    real_t* p2;
    real_t* p3;
    real_t* p4;
    real_t* p5;
    real_t* p6;
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

// TNN <-- ThermoNuclear Network

// Contain thermonuclear network data, intended to be passed between compute
// kernels. (type real_t)
struct tnn_f {
    real_t t9;
    real_t rho;
    real_t t_max;
    real_t dt_init;
};

// Contain thermonuclear network data, intended to be passed between compute
// kernels. (type real_t*)
struct tnn_fptr {
    real_t* aa;
    real_t* x;
    real_t* y;
    real_t* mass_excess;
};

// Contain thermonuclear network data, intended to be passed between compute
// kernels. (type int*)
struct tnn_iptr {
    int* z;
    int* n;
};

// Contain thermonuclear network info data.
struct tnn_info {
    int number_species;
    real_t* part_func_temp;
    char** iso_label;
    real_t** part_func;
};

// Contain one pointer for each thermonuclear network data struct.
struct tnn {
    struct tnn_info* info;
    struct tnn_f* f;
    struct tnn_fptr* fptr;
    struct tnn_iptr* iptr;
};

// Contain parameters specific to thermonuclear isotope evolution computation.
struct problem_parameters {
    int f_plus_total;
    int f_minus_total;
    real_t* f_plus;        // Dynamically-allocated 1D array for non-zero F+ (Dim
                          // f_plus_total)
    real_t* f_minus;       // Dynamically-allocated 1D array for non-zero F- (Dim
                          // f_minus_total)
    real_t* f_plus_factor; // Dynamically allocated 1D array for number species
                          // factor in F+ terms
    real_t* f_minus_factor; // Dynamically allocated 1D array for number species
                           // factor in F- terms
    real_t* f_plus_sum;     // Sum of F+ for each isotope
    real_t* f_minus_sum;    // Sum of F- for each isotope
    real_t* prefactor;
    int* f_plus_max;       // Upper index for each isotope in the f_plus array
    int* f_plus_min;       // Lower index for each isotope in the f_plus array
    int* f_minus_max;      // Upper index for each isotope in the f_minus array
    int* f_minus_min;      // Lower index for each isotope in the f_minus array
    int* f_plus_isotope_cut;  // Upper index for each isotope in Fplus (Dim
                              // number_species)
    int* f_minus_isotope_cut; // Upper index for each isotope in Fminus (Dim
                              // number_species)
    int* f_plus_number;       // Number of finite F+ components for each isotope
                              // (Dim number_species)
    int* f_minus_number;      // Number of finite F- components for each isotope
                              // (Dim number_species)
    int* f_plus_isotope_idx;  // Array containing the isotope index for each F+
                              // (Dim f_plus_total)
    int* f_minus_isotope_idx; // Array containing the isotope index for each F-
                              // (Dim f_minus_total)
    int* f_plus_map;          // Index mapper for Fplus (Dim f_plus_total)
    int* f_minus_map;         // Index mapper for Fminus (Dim f_minus_total)
    int** reaction_mask;
};

// Free all memory of source rate_library struct.
int rate_library_destroy(struct rate_library** src);

// Free all memory of source network struct.
int network_destroy(struct tnn** src);

// Free all memory of source problem_parameters struct.
int problem_parameters_destroy(struct problem_parameters** src, int species);

// Print network data.
int network_print(const struct tnn* inter_fptr);

// Print rate_library struct data.
int rate_library_print(const struct rate_library* rates,
                       const struct tnn* network);

// Return pointer to clone of source rate_library.
struct rate_library* rate_library_clone(const struct rate_library* src);

// Return pointer to clone of source tnn.
struct tnn* tnn_clone( const struct tnn* src);

#endif
