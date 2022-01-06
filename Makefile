CC = gcc

CFLAGS = -g -std=c99 -Wall

LDFLAGS = -lcmocka

.PHONY: clean all check

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

#For all compiling, building and executing
all: DNA check 

#For only building and testing interface
build: DNA DNA_bin

#Clean compilation files & output result
clean :
	@rm -rf __pycache__ build .pytest_cache output *.so
	@rm -f *.o test_gene test_gene_bin

#For only executing tests
check: run_test_gene test_DNA run_test_gene_bin test_DNA_bin

#For only running the non-binary program
run:
	sudo python3 setup.py install
	python3 main.py

#For only running the binary program
run_bin:
	sudo python3 setup_bin.py install
	python3 main_bin.py


# Naïve library
test_gene.o: gene.c

test_gene: test_gene.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

run_test_gene: test_gene
	./test_gene &

DNA : 
	python3 setup.py build
	cp build/lib*/*.so .

test_DNA : 
	python3 -m pytest -s test_DNA.py


# Binary optimized library
test_gene_bin.o: gene_bin.c

test_gene_bin: test_gene_bin.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

run_test_gene_bin: test_gene_bin
	./test_gene_bin &

DNA_bin : 
	python3 setup_bin.py build
	cp build/lib*/*.so .

test_DNA_bin : 
	python3 -m pytest -s test_DNA_bin.py
