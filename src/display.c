#include "display.h"

#include <stdio.h>
#include <stdlib.h>

int print_real_t_3d(real_t*** grid, int i_, int j_, int k_) {
    for (int i = 0; i < i_; i++) {
        for (int j = 0; j < j_; j++) {
            for (int k = 0; k < k_; k++) {
                printf(" %12.8f", grid[k][j][i]);
            }
        }
    }
    printf("\n");
    return EXIT_SUCCESS;
}

int fprint_real_t_3d(FILE* file, real_t*** grid, int i_, int j_, int k_) {
    for (int k = 0; k < k_; k++) {
        for (int j = 0; j < j_; j++) {
            for (int i = 0; i < i_; i++) {
                if (i + j + k == 0) {
                    fprintf(file, "%.8f", grid[k][j][i]);
                } else {
                    fprintf(file, " %.8f", grid[k][j][i]);
                }
            }
        }
    }
    fprintf(file, "\n");
    return EXIT_SUCCESS;
}
