#include "display.h"

#include <stdlib.h>
#include <stdio.h>

int print_float_2d(float** grid, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf(" %12.8f", grid[i][j]);
        }
    }
    printf("\n");
    return EXIT_SUCCESS;
}
