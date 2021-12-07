all: test

.PHONY: clean

run: test
	./main

test: gene.c main.c
	gcc -Wall -Wextra gene.c main.c -o  main

clean:
	rm -f *.o main