CC = gcc

CFLAGS = -g -std=c99

LDFLAGS = -lcmocka

.PHONY: clean all check

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

#For all compiling, building and executing
all: main DNA check

#For only building and testing interface
build: DNA test_DNA

main: main.o gene.o
	$(CC) $(CFLAGS) -o $@ $^

gene.o: gene.h

check: test_gene test_DNA test_gene_bin

test_gene: test_gene.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_gene.o: gene.c

test_gene_bin: test_gene_bin.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_gene_bin.o: gene_bin.c

DNA : 
	python3 setup.py build
	cp build/lib*/*.so .

test_DNA : 
	python3 -m pytest -s

clean :
	@rm -rf __pycache__ build .pytest_cache *.so
	@rm -f *.o main test_gene test_gene_bin

