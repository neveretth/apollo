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
	
APOLLO_ROCM_RUN_CMD = $(BUILD_DIR)/apollo-rocm \
	-C config/config.toml \
	-S simulation/profile.toml \
	
BUILD_DIR = $(PWD)/build
SOURCE_DIR = $(PWD)/src

default-target: apollo

clean:
	rm -rf $(BUILD_DIR)/*

apollo: builddir source kernel driver types parser toml
	$(CC) $(BUILD_DIR)/src/*.o -o $(BUILD_DIR)/apollo $(CLIBS)
	
apollo-rocm: builddir rocm-defs rocm source kernel-rocm driver types parser toml
	$(HIPCC) $(BUILD_DIR)/src/*.o -o $(BUILD_DIR)/apollo-rocm $(CLIBS)

rocm-defs: 
	$(eval CFLAGS += -D__MP_ROCM)
	
rebuild: clean apollo

builddir:
	mkdir -p $(BUILD_DIR)/src
	mkdir -p $(BUILD_DIR)/src/rocm

include Makefile.test

include src/Makefile
