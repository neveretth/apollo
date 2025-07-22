#ifndef __NUMEFFECT_H
#define __NUMEFFECT_H

#include "defs.h"

#include <stdlib.h>

// Apply randomness of one less order of magnitude than the given data.
int effect_rand(real_t*** data, int i_, int j_, int k_);

// Apply radial of one less order of magnitude than the given data.
int effect_radial(real_t*** data, int i_, int j_, int k_);

// Apply gradient of one less order of magnitude than the given data.
int effect_gradient(real_t*** data, int i_, int j_, int k_);

#endif
