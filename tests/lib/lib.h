#ifndef __UNIT_TEST_LIB_H
#define __UNIT_TEST_LIB_H

#define SCALAR 1

// Primary unit test driver. Returns EXIT_SUCCESS(FAILURE).
int unit_test_driver(int kernel_driver(void**, int*),
                     int kernel_data_init(void**, int*),
                     int kernel_data_check(void**, int*, void**, int*),
                     int num_params);

#endif
