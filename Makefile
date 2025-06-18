CC = /bin/gcc
HIPCC = /opt/rocm/bin/hipcc
CFLAGS = -std=c99 -O3 -Wall -Wextra --pedantic
HIPFLAGS = -x c
CLIBS = -lm
SRC = $(wildcard src/*.c)
HIP_SRC = $(wildcard src/*.hip)
HEADS = $(wildcard src/*.h)

default-target: apollo

clean:
	rm -rf build/*

# Link binary with HIPCC
apollo: builddir $(SRC) $(HIP_SRC)
	$(HIPCC) build/src/*.o -o build/apollo

builddir:
	mkdir -p build/src

# Compile C source
$(SRC):
	$(CC) -o build/$@.o -c $@ $(CFLAGS) $(CLIBS)
	
# Compile HIP source
$(HIP_SRC):
	$(HIPCC) -o build/$@.o -c $@

