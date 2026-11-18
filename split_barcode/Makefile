CC      ?= gcc
CFLAGS  ?= -O3 -march=native -std=gnu99 -Wall -Wextra
LIBS     = -lz -lpthread

all: split_barcode

split_barcode: split_barcode.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f split_barcode

.PHONY: all clean
