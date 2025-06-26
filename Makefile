CC = /bin/gcc
HIPCC = /opt/rocm/bin/hipcc
CFLAGS = -std=c99 -O3 -Wall -Wextra --pedantic -g3
HIPFLAGS = -x c
CLIBS = -lm -lhdf5
SRC = $(wildcard src/*.c)
HIP_SRC = $(wildcard src/*.hip)
HEADS = $(wildcard src/*.h)

PROFILEFLAGS = CPUPROFILE=/tmp/prof.out

default-target: apollo

clean:
	rm -rf build/*

test: apollo
	/bin/time \
	-f "\n| Time: %e | CPU: %P | User: %U | Kernel: %S | MinPF: %R | MajPF: %F |" \
	./build/apollo \
	--hydro-debug \
	--neutrino-file data/FENNData40M186.h5 \
	--rate-library data/rateLibrary_alpha.data \
	--network data/CUDAnet_alpha.inp \
	# --verbose \
	# --thermo-debug \
	# --neutrino-debug \
	
test-mem: apollo
	valgrind \
	./build/apollo \
	--hydro-debug \
	--neutrino-file data/FENNData40M186.h5 \
	--rate-library data/rateLibrary_alpha.data \
	--network data/CUDAnet_alpha.inp \
	# --verbose \
	# --thermo-debug \
	# --neutrino-debug \

test-prof: apollo
	$(HIPCC) build/src/*.o -o build/apollo-prof $(CLIBS) -lprofiler
	
	$(PROFILEFLAGS) ./build/apollo-prof \
	--hydro-debug \
	--neutrino-file data/FENNData40M186.h5 \
	--rate-library data/rateLibrary_alpha.data \
	--network data/CUDAnet_alpha.inp
	
	pprof ./build/apollo-prof /tmp/prof.out

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
