CC = /bin/gcc
HIPCC = /opt/rocm/bin/hipcc
CFLAGS = -std=c99 -O3 -Wall -Wextra --pedantic -g3
HIPFLAGS = -x c
CLIBS = -lm -lhdf5

APOLLO_RUN_CMD = $(BUILD_DIR)/apollo \
	--neutrino-file data/FENNData40M186.h5 \
	--rate-library data/rateLibrary_alpha.data \
	--network data/CUDAnet_alpha.inp \
	--full \
	# --neutrino-debug \
	# --hydro-debug \
	# --verbose \
	# --thermo-debug \
	
BUILD_DIR = $(PWD)/build
SOURCE_DIR = $(PWD)/src

default-target: apollo

clean:
	rm -rf $(BUILD_DIR)/*

apollo: builddir source kernel driver types
	$(HIPCC) $(BUILD_DIR)/src/*.o -o $(BUILD_DIR)/apollo $(CLIBS)
	
builddir:
	mkdir -p $(BUILD_DIR)/src

include Makefile.test

include src/Makefile
