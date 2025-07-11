CC = /bin/gcc
HIPCC = /opt/rocm/bin/hipcc
CFLAGS = -std=c99 -O3 -Wall -Wextra --pedantic -g3
HIPFLAGS = -x c
CLIBS = -lm -lhdf5

# Comment out TIME_CMD if you don't want to time Apollo
TIME_CMD = /bin/time -f "\n| Time: %e | CPU: %P | User: %U | Kernel: %S | MinPF: %R | MajPF: %F |"

APOLLO_RUN_CMD = $(BUILD_DIR)/apollo \
	--rate-library data/ratelibrary-alpha.aad \
	--network data/network-alpha.aad \
	--neutrino-file data/FENNData40M186.h5 \
	--hydro-debug \
	# --thermo-debug \
	# --neutrino-debug \
	# --full \
	# --verbose \
	
APOLLO_ROCM_RUN_CMD = $(BUILD_DIR)/apollo-rocm \
	--rate-library data/ratelibrary-alpha.aad \
	--network data/network-alpha.aad \
	--neutrino-file data/FENNData40M186.h5 \
	--thermo-debug \
	--rocm-accel \
	# --neutrino-debug \
	# --full \
	# --hydro-debug \
	# --verbose \
	
BUILD_DIR = $(PWD)/build
SOURCE_DIR = $(PWD)/src

default-target: apollo apollo-rocm

clean:
	rm -rf $(BUILD_DIR)/*

apollo: builddir source kernel driver types parser
	$(CC) $(BUILD_DIR)/src/*.o -o $(BUILD_DIR)/apollo $(CLIBS)
	
apollo-rocm: builddir source-rocm kernel-rocm driver-rocm types parser
	@# shitty hack to remove duplicate symbol.
	rm -rf $(BUILD_DIR)/src/kernel-driver.c.o
	$(HIPCC) $(BUILD_DIR)/src/*.o $(BUILD_DIR)/src/rocm/*.o -o $(BUILD_DIR)/apollo-rocm $(CLIBS)
	
builddir:
	mkdir -p $(BUILD_DIR)/src
	mkdir -p $(BUILD_DIR)/src/rocm

include Makefile.test

include src/Makefile
