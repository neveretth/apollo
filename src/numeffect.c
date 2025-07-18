#include "numeffect.h"

#include <math.h>
#include <stdio.h>

int effect_rand(float*** data, int i_, int j_, int k_) {
    for (int i = 0; i < i_; i++) {
        for (int j = 0; j < j_; j++) {
            for (int k = 0; k < k_; k++) {
                data[i][j][k] +=
                    (float)rand() / (10 * (float)(RAND_MAX / data[i][j][k]));
            }
        }
    }
    return EXIT_SUCCESS;
}

int effect_radial(float*** data, int i_, int j_, int k_) {
    float rad = sqrt((i_ * i_) + (j_ * j_) + (k_ * k_));
    for (int i = 0; i < i_; i++) {
        for (int j = 0; j < j_; j++) {
            for (int k = 0; k < k_; k++) {
                float len = sqrt((i * i) + (j * j) + (k * k)) / rad;
                len = 1 - len;
                // Should this be here, before the len = 1 - len, or never.
                len *= len;
                data[i][j][k] += (data[i][j][k] / 10) * len;
            }
        }
    }
    return EXIT_SUCCESS;
}

int effect_gradient(float*** data, int i_, int j_, int k_) {
    float bound = i_;
    for (int i = 0; i < i_; i++) {
        for (int j = 0; j < j_; j++) {
            for (int k = 0; k < k_; k++) {
                data[i][j][k] += (data[i][j][k] / 10) * (i / bound);
            }
        }
    }
    return EXIT_SUCCESS;
}
