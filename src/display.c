#include "display.h"

#include <stdio.h>
#include <stdlib.h>

int print_float_2d(float** grid, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf(" %12.8f", grid[i][j]);
        }
    }
    printf("\n");
    return EXIT_SUCCESS;
}

int print_float_3d(float*** grid, int i_, int j_, int k_) {
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
