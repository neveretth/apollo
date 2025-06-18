CC = /bin/gcc
CFLAGS = -std=c99 -O3 -Wall -Wextra --pedantic
CLIBS = -lm
SRC = $(wildcard src/*.c)
HEADS = $(wildcard src/*.h)

apollo: $(SRC)
	mkdir -p build
	$(CC) -o build/$@ $^ $(CFLAGS) $(CLIBS)
