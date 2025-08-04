#include "display.h"

#include <stdio.h>
#include <stdlib.h>

int print_real_t_1d(real_t* grid, int len) {
    for (int i = 0; i < len; i++) {
        printf(" %12.8f", grid[i]);
    }
    printf("\n");
    return EXIT_SUCCESS;
}

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

int fprint_real_t_1d(FILE* file, real_t* grid, int len) {
    for (int i = 0; i < len; i++) {
        if (i == 0) {
            fprintf(file, "%12.8f", grid[i]);
        } else {
            fprintf(file, " %12.8f", grid[i]);
        }
    }
    fprintf(file, "\n");
    return EXIT_SUCCESS;
}

int fprint_real_t_2d(FILE* file, real_t** grid, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (rows + cols == 0) {
                fprintf(file, "%12.8f", grid[i][j]);
            } else {
                fprintf(file, " %12.8f", grid[i][j]);
            }
        }
    }
    fprintf(file, "\n");
    return EXIT_SUCCESS;
}

int fprint_real_t_3d(FILE* file, real_t*** grid, int i_, int j_, int k_) {
    for (int i = 0; i < i_; i++) {
        for (int j = 0; j < j_; j++) {
            for (int k = 0; k < k_; k++) {
                if (i + j + k == 0) {
                    fprintf(file, "%.8f", grid[i][j][k]);
                } else {
                    fprintf(file, " %.8f", grid[i][j][k]);
                }
            }
        }
    }
    fprintf(file, "\n");
    return EXIT_SUCCESS;
}
