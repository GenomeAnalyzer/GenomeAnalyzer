CC = gcc

CFLAGS = -g -std=c99

LDFLAGS = -lcmocka

.PHONY: clean all

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

all: main

main: main.o gene.o
	$(CC) $(CFLAGS) -o $@ $^

gene.o: gene.h

check: test_gene

test_gene: test_gene.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_gene.o: gene.c 

clean:
	rm -f *.o main test_gene