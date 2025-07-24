#!/bin/bash

# The first arg is the build directory.

printf "\n\n[test | HYDRO]\n\n"
./tests/hydro-test/test.sh $1

printf "\n\n[test | THERMO]\n\n"
./tests/thermo-test/test.sh $1

printf "\n\n[test | NEUTRINO]\n\n"
./tests/neutrino-test/test.sh $1
