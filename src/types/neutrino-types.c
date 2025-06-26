#include "neutrino-types.h"
#include "../types.h"

#include <stdlib.h>

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
