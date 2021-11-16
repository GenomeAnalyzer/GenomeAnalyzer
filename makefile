all: test

.PHONY: clean

run: test
	./main

test: gene.c main.c
	gcc -Wall -Wextra *.c -o  main

clean:
	rm -f *.o main