CC ?= gcc
CFLAGS ?= -O3 -march=native -fopenmp -Wall -Wno-unused-result
build/axi2phase: axi2phase.c
	mkdir -p build
	$(CC) $(CFLAGS) -o $@ $< -lm
clean:
	rm -rf build
