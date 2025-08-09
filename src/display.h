#ifndef __DISPLAY_H
#define __DISPLAY_H

#include "defs.h"

#include <stdio.h>

// Print 3D grid of floating point numbers to stdout.
int print_real_t_3d(real_t*** grid, int i, int j, int k);

// Print 3D grid of floating point numbers to file.
int fprint_real_t_3d(FILE* file, real_t*** grid, int i, int j, int k);

#endif
