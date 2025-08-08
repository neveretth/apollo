#include "lib.h"

#include <errno.h>
#include <fcntl.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define VALIDATE_STRING "validate"
#define VALIDATE_STRLEN 8

static int dump_data(char* filepath, void** param, int* size, int length) {
    bool fail = false;

    int dumpfd = open(filepath, O_RDWR);
    if (errno != 0) {
        perror(filepath);
        fail = true;
        goto exit;
    }
    // Start by writing a validation string. (no null terminus)
    char validate[VALIDATE_STRLEN] = VALIDATE_STRING;
    int validate_strlen = VALIDATE_STRLEN;
    int err;
    err = write(dumpfd, validate, VALIDATE_STRLEN * sizeof(char));
    if (err == -1) {
        fail = true;
        goto exit;
    }
    // For the rest of the time I'm going to assume this call doesn't
    // fail to simplify the code.

    // Write the len of the size ptr (of type int, size i32)
    // we really should use <stdint.h>...
    err = write(dumpfd, &length, sizeof(int));

    // Write the sizes
    err = write(dumpfd, size, length * sizeof(int));

    // Write the data
    for (int i = 0; i < length; i++) {
        err = write(dumpfd, param[i], size[i]);
    }

    // XOR the validate string for simple corruption checking.
    int byte = 0;
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < size[i]; j++) {
            // Validate by XORing each byte of the parameters against
            // sequential bytes of the validate string.
            validate[byte % validate_strlen] ^= ((char*)param[i])[j];
            byte++;
        }
    }
    // Finally, write the xor-ed validate string.
    err = write(dumpfd, validate, VALIDATE_STRLEN * sizeof(char));

    // Check for any errors during the process.
    if (err == -1) {
        perror("");
        fail = true;
        goto exit;
    }

exit:
    close(dumpfd); // This may cause issues if it did not open.
    if (fail) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

// This function is fundamentally similar to dump_data.
static int read_data(char* filepath, void** param, int* size) {
    bool fail = false;

    int readfd = open(filepath, O_RDONLY);
    if (errno != 0) {
        perror(filepath);
        fail = true;
        goto exit;
    }

    char validate[VALIDATE_STRLEN];
    int validate_strlen = VALIDATE_STRLEN;
    int err;
    err = read(readfd, validate, VALIDATE_STRLEN * sizeof(char));
    if (err == -1) {
        fail = true;
        goto exit;
    }

    if (strncmp(validate, VALIDATE_STRING, VALIDATE_STRLEN) != 0) {
        printf("[APOLLO UNIT TEST] Data read validation init failed at %s\n",
               filepath);
        fail = true;
        goto exit;
    }

    int length;

    err = read(readfd, &length, sizeof(int));

    err = read(readfd, size, length * sizeof(int));

    for (int i = 0; i < length; i++) {
        err = read(readfd, param[i], size[i]);
    }

    int byte = 0;
    for (int i = 0; i < length; i++) {
        for (int j = 0; j < size[i]; j++) {
            validate[byte % validate_strlen] ^= ((char*)param[i])[j];
            byte++;
        }
    }
    char tmp[VALIDATE_STRLEN];
    err = read(readfd, tmp, VALIDATE_STRLEN * sizeof(char));

    if (strncmp(validate, tmp, VALIDATE_STRLEN) != 0) {
        printf("[APOLLO UNIT TEST] Data read validation check failed at %s\n",
               filepath);
        fail = true;
        goto exit;
    }

    if (err == -1) {
        perror("");
        fail = true;
        goto exit;
    }

exit:
    close(readfd);
    if (fail) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int unit_test_driver(int kernel_driver(void**, int*),
                     int kernel_data_init(void**, int*),
                     int kernel_data_check(void**, int*, void**, int*),
                     int num_params) {

    // Track possible error occurance.
    bool fail = false;

    // Initialize the data to be read. While this does assign values, they will
    // be overwritten. The only way this could give inaccurate results is if the
    // kernel does not do anything. (I am not against improving this though...)
    void** param_read = malloc(num_params * sizeof(void*));
    int* size_read = malloc(num_params * sizeof(int));
    if (kernel_data_init(param_read, size_read) == EXIT_FAILURE) {
        fail = true;
        goto exit;
    }

    // Temporary storage of snapshot data.
    char filepath[] = "tests/sample/snapshot.data";

    if (read_data(filepath, param_read, size_read) == EXIT_FAILURE) {
        fail = true;
        goto exit;
    }
    
    // Initialize the data for the kernel test.
    void** param = malloc(num_params * sizeof(void*));
    int* size = malloc(num_params * sizeof(int));
    if (kernel_data_init(param, size) == EXIT_FAILURE) {
        fail = true;
        goto exit;
    }

    // For future versions, we may read data here as well before
    // the kernel. That is trivial to implement.

    // Run the kernel with the initialized params.
    if (kernel_driver(param, size) == EXIT_FAILURE) {
        return EXIT_FAILURE;
    }

    // Some sort of comparison between the two sets of data...
    if (kernel_data_check(param_read, size_read, param, size) == EXIT_FAILURE) {
        fail = true;
        goto exit;
    }

    // This exists for generating data, but we do not want to overwrite
    // it on normal runs. Eventually this will be made into an option
    // given, perhaps: unit-test --dump

    // dump_data(filepath, param, size, num_params);

exit:
    // Free the allocated memory.
    for (int i = 0; i < num_params; i++) {
        free(param[i]);
        free(param_read[i]);
    }
    free(size);
    free(param);
    free(size_read);
    free(param_read);

    if (fail) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
