CC=gcc -O2 -g -std=gnu99

.PHONY: all
.SUFFIXES: .c .o

APPS=bem_convert bem_generate

all: $(APPS)

.c.o:
	$(CC) -c $< -o $@

bem_convert: bem_convert.o bem_file.o bem_aux.o
	$(CC) $^ -o $@

bem_generate: bem_generate.o bem_file.o bem_aux.o
	$(CC) $^ -o $@

bem_test: bem_file.c
	$(CC) -DTEST bem_file.c -o $@

clean:
	rm -rf *.o $(APPS)
