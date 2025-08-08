#include "lib.h"

#include <stdio.h>
#include <stdlib.h>

// Only included because we are hardcoding for the sample kernel
#include "../sample/link.h"
#include "../sample/logic.h"

int main(int argc, char** argv) {

    printf("[APOLLO UNIT TEST] Beginning unit tests...\n");

    // For the time being we will hardcode the behavior.
    // Passing functions like this allows us to keep the test logic
    // itself entirely kernel-agnostic. Only the drivers and data init
    // need to know anything about the kernel.

    printf("[APOLLO UNIT TEST] Testing good kernel...\n");
    if (unit_test_driver(sample_kernel_driver, sample_kernel_data_init,
                         sample_kernel_data_check,
                         SAMPLE_KERNEL_NUM_PARAMS) == EXIT_FAILURE) {
        printf("[APOLLO UNIT TEST] Unit test failed.\n");
        return EXIT_FAILURE;
    }
    printf("[APOLLO UNIT TEST] Test passed :D\n");

    // To see it pass, comment out the code between the dividers.
    //////////// \\\\\\\\\\\\;
    printf("[APOLLO UNIT TEST] Testing bad kernel...\n");
    if (unit_test_driver(sample_kernel_bad_driver, sample_kernel_data_init,
                         sample_kernel_data_check,
                         SAMPLE_KERNEL_NUM_PARAMS) == EXIT_FAILURE) {
        printf("[APOLLO UNIT TEST] Unit test failed.\n");
        return EXIT_FAILURE;
    }
    printf("[APOLLO UNIT TEST] Test passed :D\n");
    // \\\\\\\\\\\\ ////////////;

    printf("[APOLLO UNIT TEST] All tests passed ~~ Exiting.\n");
    return EXIT_SUCCESS;
}
