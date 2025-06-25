CC = /bin/gcc
HIPCC = /opt/rocm/bin/hipcc
CFLAGS = -std=c99 -O3 -Wall -Wextra --pedantic
HIPFLAGS = -x c
CLIBS = -lm -lhdf5
SRC = $(wildcard src/*.c)
HIP_SRC = $(wildcard src/*.hip)
HEADS = $(wildcard src/*.h)

default-target: apollo

clean:
	rm -rf build/*

test: apollo
	./build/apollo \
	--neutrino-file data/FENNData40M186.h5 \
	--neutrino-debug \
	# --rate-library data/rateLibrary_alpha.data \
	# --network data/CUDAnet_alpha.inp \
	# --thermo-debug \

# Link binary with HIPCC
apollo: builddir $(SRC) $(HIP_SRC)
	$(HIPCC) build/src/*.o -o build/apollo $(CLIBS)

builddir:
	mkdir -p build/src

# Compile C source
$(SRC):
	$(CC) -o build/$@.o -c $@ $(CFLAGS) $(CLIBS)
	
# Compile HIP source
$(HIP_SRC):
	$(HIPCC) -o build/$@.o -c $@
