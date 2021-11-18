CC = gcc

CFLAGS = -g -std=c99

LDFLAGS = -lcmocka

.PHONY: clean all check

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

all: main DNA check

main: main.o gene.o
	$(CC) $(CFLAGS) -o $@ $^

gene.o: gene.h

check: test_gene test_DNA

test_gene: test_gene.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_gene.o: gene.c 

DNA : 
	python3 setup.py build
	cp build/lib*/*.so .

test_DNA : 
	python3 -m pytest

clean :
	@rm -rf __pycache__ build .pytest_cache *.so
	rm -f *.o main test_gene

