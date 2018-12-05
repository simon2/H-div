CC=gcc -O2 -g -std=gnu99

.PHONY: all
.SUFFIXES: .c .o

all: bem_convert

.c.o:
	$(CC) -c $< -o $@

bem_convert: bem_convert.o bem_file.o
	$(CC) $^ -o $@

bem_test: bem_file.c
	$(CC) -DTEST bem_file.c -o $@
