CC = gcc

CFLAGS = -g -std=c99

.PHONY: clean all check

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

all: run run_bin

run: main.c ../gene.c
	$(CC) $(CFLAGS) -o $@ $^
	
run_bin: main_bin.c ../gene_bin.c
	$(CC) $(CFLAGS) -o $@ $^

clean :
	@rm -f *.o run run_bin