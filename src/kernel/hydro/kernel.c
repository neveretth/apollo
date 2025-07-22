#include "kernel.h"

#include <math.h>
#include <stdlib.h>

int hydro_integration_kernel(float*** temp, float*** density, float volume,
                                float h, float dt, float t_end, int* dim) {
    // This really should be passed, not created each time.
    float*** delta_temp = malloc(dim[0] * sizeof(float*));
    for (int i = 0; i < dim[0]; i++) {
        delta_temp[i] = malloc(dim[1] * sizeof(float*));
        for (int j = 0; j < dim[1]; j++) {
            delta_temp[i][j] = malloc(dim[2] * sizeof(float));
        }
    }

    // For now just assume it's all even. Volume is for each cell.
    float area = volume / dim[2]; // Placeholder area of interaction.

    float c = 1e8;  // Placeholder contribution value.
    float ntd = 0;  // Nuclear burning temp diff (assume negligible for now)

    float t = 0;

    // Does one integration (single step if t_end = 0)
    do {
        for (int i = 0; i < dim[0]; i++) {
            for (int j = 0; j < dim[1]; j++) {
                for (int k = 0; k < dim[2]; k++) {
                    float pdvh = 1 / (density[i][j][k] * volume * c);
                    float temp_diff = -(6 * temp[i][j][k]);

                    // As with the linear integration, towards idx=0 is positive
                    temp_diff += temp[i - (i > 0)][j][k];
                    temp_diff += temp[i + (i < dim[0] - 1)][j][k];
                    temp_diff += temp[i][j - (j > 0)][k];
                    temp_diff += temp[i][j + (j < dim[1] - 1)][k];
                    temp_diff += temp[i][j][k - (k > 0)];
                    temp_diff += temp[i][j][k + (k < dim[2] - 1)];

                    delta_temp[i][j][k] =
                        dt * (pdvh * (ntd + (h * area * temp_diff)));
                }
            }
        }
        for (int i = 0; i < dim[0]; i++) {
            for (int j = 0; j < dim[1]; j++) {
                for (int k = 0; k < dim[2]; k++) {
                    temp[i][j][k] += delta_temp[i][j][k];
                }
            }
        }
        t += dt;
    } while (t < t_end);

    for (int i = 0; i < dim[0]; i++) {
        for (int j = 0; j < dim[1]; j++) {
            free(delta_temp[i][j]);
        }
        free(delta_temp[i]);
    }
    free(delta_temp);

    return EXIT_SUCCESS;
}

int hydro_data_preprocess() {
    return EXIT_SUCCESS;
}
