CC = /bin/gcc
HIPCC = /opt/rocm/bin/hipcc
CFLAGS = -std=c99 -O3 -Wall -Wextra --pedantic -g3
HIPFLAGS = -x c
CLIBS = -lm -lhdf5

# Comment out TIME_CMD if you don't want to time Apollo
TIME_CMD = /bin/time -f "\n| Time: %e | CPU: %P | User: %U | Kernel: %S | MinPF: %R | MajPF: %F |"

APOLLO_RUN_CMD = $(BUILD_DIR)/apollo \
	--rate-library data/ratelibrary-365.aad \
	--network data/network-365.aad \
	--thermo-debug \
	# --neutrino-file data/FENNData40M186.h5 \
	# --full \
	# --neutrino-debug \
	# --hydro-debug \
	# --verbose \
	
BUILD_DIR = $(PWD)/build
SOURCE_DIR = $(PWD)/src

default-target: apollo

clean:
	rm -rf $(BUILD_DIR)/*

apollo: builddir source kernel driver types parser
	$(HIPCC) $(BUILD_DIR)/src/*.o -o $(BUILD_DIR)/apollo $(CLIBS)
	
builddir:
	mkdir -p $(BUILD_DIR)/src

include Makefile.test

include src/Makefile
