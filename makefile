all: test

.PHONY: clean

test: gene.c main.c
	gcc -Wall -Wextra *.c -o  main

clean:
	rm -f *.o main