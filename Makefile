CC = /bin/gcc
HIPCC = /opt/rocm/bin/hipcc
CFLAGS = -std=c99 -O3 -Wall -Wextra --pedantic -g3
HIPFLAGS = -x c
CLIBS = -lm -lhdf5 -lc

# Comment out TIME_CMD if you don't want to time Apollo
TIME_CMD = /bin/time -f "\n| Time: %e | CPU: %P | User: %U | Kernel: %S | MinPF: %R | MajPF: %F |"

APOLLO_RUN_CMD = $(BUILD_DIR)/apollo \
	-C config/config.toml \
	-S simulation/profile.toml \
	# --hydro-debug \
	# --rate-library data/ratelibrary-alpha.aad \
	# --network data/network-alpha.aad \
	# --neutrino-file data/FENNData40M186.h5 \
	# --neutrino-debug \
	# --full \
	# --thermo-debug \
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

default-target: apollo

clean:
	rm -rf $(BUILD_DIR)/*

apollo: builddir source kernel driver types parser toml
	$(CC) $(BUILD_DIR)/src/*.o -o $(BUILD_DIR)/apollo $(CLIBS)
	
apollo-rocm: builddir source-rocm kernel-rocm driver-rocm types parser toml
	@# shitty hack to remove duplicate symbol.
	rm -rf $(BUILD_DIR)/src/kernel-driver.c.o
	$(HIPCC) $(BUILD_DIR)/src/*.o $(BUILD_DIR)/src/rocm/*.o -o $(BUILD_DIR)/apollo-rocm $(CLIBS)
	
rebuild: clean apollo

builddir:
	mkdir -p $(BUILD_DIR)/src
	mkdir -p $(BUILD_DIR)/src/rocm

include Makefile.test

include src/Makefile
