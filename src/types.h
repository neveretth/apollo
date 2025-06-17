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
    float* rate;
    float* flux;
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

// Contain network data.
struct network {
    int* z;
    int* n;
    float* aa;
    float* x;
    float* y;
    float* mass_excess;
    float* part_func_temp;
    float** part_func;
    char** iso_label;
    int number_species;
};

// Contain problem parameters,
// TODO: rename this.
struct problem_parameters {
    int f_plus_total;
    int f_minus_total;
    float* f_plus;        // Dynamically-allocated 1D array for non-zero F+ (Dim
                          // f_plus_total)
    float* f_minus;       // Dynamically-allocated 1D array for non-zero F- (Dim
                          // f_minus_total)
    float* f_plus_factor; // Dynamically allocated 1D array for number species
                          // factor in F+ terms
    float* f_minus_factor; // Dynamically allocated 1D array for number species
                           // factor in F- terms
    float* f_plus_sum;     // Sum of F+ for each isotope
    float* f_minus_sum;    // Sum of F- for each isotope
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
int rate_library_destroy(struct rate_library* src);

// Free all memory of source network struct.
int network_destroy(struct network* src);

// Free all memory of source problem_parameters struct.
int problem_parameters_destroy(struct problem_parameters* src, int species);

// Print network struct data.
int network_print(const struct network* network);

// Print rate_library struct data.
int rate_library_print(const struct rate_library* rates,
                       const struct network* network);

// Print abundances present in network.
int print_abundances(const struct network* network);

// Print result data.
int print_results(const struct rate_library* rates,
                  const struct network* network,
                  const struct problem_parameters* params);

#endif
