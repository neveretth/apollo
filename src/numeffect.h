#ifndef __NUMEFFECT_H
#define __NUMEFFECT_H

#include <stdlib.h>

// Apply randomness of one less order of magnitude than the given data.
int effect_rand(float*** data, int i_, int j_, int k_);

// Apply radial of one less order of magnitude than the given data.
int effect_radial(float*** data, int i_, int j_, int k_);

// Apply randomness of one less order of magnitude than the given data.
int effect_gradient(float*** data, int i_, int j_, int k_);

#endif
