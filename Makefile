# Define all vars that are not already defined.
CC ?= /bin/gcc
HIPCC ?= /opt/rocm/bin/hipcc
# GCCOPTS is for optimization. We can do away with
# much of strict compliance, but we _cannot_ use
# -ffast-math as it breaks too much IEEE 754 and while
# this is _probably_ fine, we choose to avoid it.
# Essentially, we assume apollo is being well behaved (note
# that it is often _not_, we just _assume) and optimize for
# that. If the end result is nan/inf... whoops, maybe comment
# out GCCOPTS and try again (or make changes to the args).
GCCOPTS ?= -O3 -fno-math-errno -ffinite-math-only -frounding-math
CFLAGS ?= -std=c99 -Wall -Wextra --pedantic -g3 -D__REAL_TYPE_FLOAT $(GCCOPTS)
HIPFLAGS ?= -x c
CLIBS ?= -lhdf5 -lhdf5_hl -lc -lm

# Custom user make variable setup.
-include Makefile.user

# Export user-overrides so sub-makefiles can see them
export CLIBS
export CFLAGS
export LDFLAGS


# Comment out TIME_CMD if you don't want to time Apollo
TIME_CMD = /bin/time -f "\n| Time: %e | CPU: %P | User: %U | Kernel: %S | MinPF: %R | MajPF: %F |"

APOLLO_RUN_CMD = $(BUILD_DIR)/apollo \
	-C config/config.toml \
	-S simulation/simulation.toml \
	
APOLLO_ROCM_RUN_CMD = $(BUILD_DIR)/apollo-rocm \
	-C config/config.toml \
	-S simulation/simulation.toml \
	
BUILD_DIR = $(PWD)/build
SOURCE_DIR = $(PWD)/src

default-target: apollo

clean:
	rm -rf $(BUILD_DIR)/src

apollo: builddir source kernel driver types parser toml
	$(CC) $(BUILD_DIR)/src/*.o -o $(BUILD_DIR)/apollo $(CLIBS)
	
apollo-rocm: rocm-defs builddir rocm source kernel-rocm driver types parser toml
	$(HIPCC) $(BUILD_DIR)/src/*.o -o $(BUILD_DIR)/apollo-rocm $(LDFLAGS) $(CLIBS)

rocm-defs: 
	$(eval CFLAGS += -D__MP_ROCM)
	
rebuild: clean apollo

builddir:
	mkdir -p $(BUILD_DIR)/src
	mkdir -p $(BUILD_DIR)/src/rocm

include Makefile.test

include src/Makefile
