#ifndef __VALIDATE_H
#define __VALIDATE_H

#include <stdio.h>

// Return status EXIT_SUCCESS(FAILURE) after checking file. Print diagnostic
// information if failure occurs.
int validate_file(FILE* file);

#endif
