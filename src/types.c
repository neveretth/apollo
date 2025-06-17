#include "types.h"

#include <stdlib.h>

// Don't bother setting to NULL becuase they're all inaccessible anyways.
int rate_library_destroy(struct rate_library* src) {
    free(src->rg_class);
    free(src->rg_member_idx);
    free(src->reaction_library_class);
    free(src->num_react_species);
    free(src->num_products);
    free(src->is_ec);
    free(src->is_reverse);
    free(src->prefactor);
    free(src->q_value);
    free(src->p0);
    free(src->p1);
    free(src->p2);
    free(src->p3);
    free(src->p4);
    free(src->p5);
    free(src->p6);
    free(src->reactant_1);
    free(src->reactant_2);
    free(src->reactant_3);
    for (int i = 0; i < src->number_reactions; i++) {
        free(src->reactant_z[i]);
        free(src->reactant_n[i]);
        free(src->reactant_idx[i]);
        free(src->product_z[i]);
        free(src->product_n[i]);
        free(src->product_idx[i]);
        free(src->reaction_label[i]);
    }
    free(src->reactant_z);
    free(src->reactant_n);
    free(src->product_z);
    free(src->reactant_idx);
    free(src->product_n);
    free(src->product_idx);
    free(src->reaction_label);
    free(src);
    return EXIT_SUCCESS;
}

int network_destroy(struct network* src) {
    free(src->z);
    free(src->n);
    free(src->aa);
    free(src->y);
    free(src->x);
    free(src->mass_excess);
    for (int i = 0; i < src->number_species; i++) {
        free(src->iso_label[i]);
        free(src->part_func[i]);
    }
    free(src->iso_label);
    free(src->part_func);
    free(src);
    return EXIT_SUCCESS;
}
