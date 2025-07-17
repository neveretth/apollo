#ifndef __DISPLAY_H
#define __DISPLAY_H

#include <stdio.h>

// Print 2d grid of floating point numbers to stdout.
int print_float_2d(float** grid, int rows, int cols);

int print_float_3d(float*** grid, int i, int j, int k);

// Print grid to file
int fprint_float_3d(FILE* file, float*** grid, int i, int j, int k);

#endif
