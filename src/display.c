#include "display.h"

#include <stdio.h>
#include <stdlib.h>

int print_real_t_2d(real_t** grid, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf(" %12.8f", grid[i][j]);
        }
    }
    printf("\n");
    return EXIT_SUCCESS;
}

int print_real_t_3d(real_t*** grid, int i_, int j_, int k_) {
    for (int i = 0; i < i_; i++) {
        for (int j = 0; j < j_; j++) {
            for (int k = 0; k < k_; k++) {
                printf(" %12.8f", grid[i][j][k]);
            }
        }
    }
    printf("\n");
    return EXIT_SUCCESS;
}

// Print grid to file
int fprint_real_t_3d(FILE* file, real_t*** grid, int i_, int j_, int k_) {
    for (int i = 0; i < i_; i++) {
        for (int j = 0; j < j_; j++) {
            for (int k = 0; k < k_; k++) {
                if (i + j + k == 0) {
                    fprintf(file, "%12.8f", grid[i][j][k]);
                } else {
                    fprintf(file, " %12.8f", grid[i][j][k]);
                }
            }
        }
    }
    fprintf(file, "\n");
    return EXIT_SUCCESS;
}
