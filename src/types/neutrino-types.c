#include "neutrino-types.h"
#include "../types.h"

#include <stdlib.h>
#include <string.h>

int neunet_destroy(struct neunet** __src) {
    struct neunet* src = *__src;
    freenptr((void*)src->fptr, 4);
    free(src->fptr);
    freenptr((void*)src->info->rate_in, src->info->num_groups);
    free(src->info->rate_in);
    freenptr((void*)src->info->rate_out, src->info->num_groups);
    free(src->info->rate_out);
    free(src->info);
    free(src->f);
    free(src);
    *__src = NULL;
    return EXIT_SUCCESS;
}

int neunet_print(struct neunet* src) {
    for (int i = 0; i < src->info->num_groups; i++) {
        printf("%f ", src->fptr->n_old[i]);
    }
    printf("\n");
    return EXIT_SUCCESS;
}

struct neunet* neunet_clone(const struct neunet* src) {
    struct neunet* dest = malloc(sizeof(struct neunet));

    dest->f = malloc(sizeof(struct neunet_f));
    memcpy(dest->f, src->f, 17 * sizeof(real_t));

    dest->fptr = malloc(sizeof(struct neunet_fptr));
    dest->fptr->dv = malloc(src->info->num_groups * sizeof(real_t));
    dest->fptr->ec = malloc(src->info->num_groups * sizeof(real_t));
    dest->fptr->n_eq = malloc(src->info->num_groups * sizeof(real_t));
    dest->fptr->n_old = malloc(src->info->num_groups * sizeof(real_t));
    memcpy(dest->fptr->dv, src->fptr->dv,
           src->info->num_groups * sizeof(real_t));
    memcpy(dest->fptr->ec, src->fptr->ec,
           src->info->num_groups * sizeof(real_t));
    memcpy(dest->fptr->n_old, src->fptr->n_old,
           src->info->num_groups * sizeof(real_t));
    memcpy(dest->fptr->n_eq, src->fptr->n_eq,
           src->info->num_groups * sizeof(real_t));

    dest->info = malloc(sizeof(struct neunet_info));
    dest->info->num_groups = src->info->num_groups;
    dest->info->rate_in = malloc(src->info->num_groups * sizeof(real_t*));
    dest->info->rate_out = malloc(src->info->num_groups * sizeof(real_t*));
    for (int i = 0; i < src->info->num_groups; i++) {
        dest->info->rate_in[i] = malloc(src->info->num_groups * sizeof(real_t));
        memcpy(dest->info->rate_in[i], src->info->rate_in[i],
               src->info->num_groups * sizeof(real_t));
        dest->info->rate_out[i] = malloc(src->info->num_groups * sizeof(real_t));
        memcpy(dest->info->rate_out[i], src->info->rate_out[i],
               src->info->num_groups * sizeof(real_t));
    }

    return dest;
}
