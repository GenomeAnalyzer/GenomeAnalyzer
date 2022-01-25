CC = gcc

CFLAGS = -g -std=c99 -Wall

LDFLAGS = -lcmocka

V0 = ./versions/naive

V1 = ./versions/bin

BIN = ./bin/

OUTPUT = ./output/

.PHONY: clean all check


ifeq (run, $(firstword $(MAKECMDGOALS)))
  runargs := $(wordlist 2, $(words $(MAKECMDGOALS)), $(MAKECMDGOALS))
  $(eval $(runargs):;@true)
endif

#For all compiling, building and executing
all: DNA check 

#For only building and testing interface
build: DNA DNA_bin

#Clean compilation files & output result
clean :
	find . -type d  -name "__pycache__" -exec rm -rv {} +
	@rm -rf  build .pytest_cache output *.so
	@rm -f $(OUTPUT)*.o $(OUTPUT)test_gene $(OUTPUT)test_gene_bin

#For only executing tests
check: run_test_gene test_DNA run_test_gene_bin test_DNA_bin

#For only running the non-binary program
run:
	sudo python3 $(V0)/python/setup.py install
	python3 $(V0)/python/main.py $(runargs)
	

#For only running the binary program
run_bin:
	sudo python3 $(V1)/python/setup_bin.py install
	python3 $(V1)/python/main_bin.py $(runargs)


# Naïve library
test_gene.o: gene

test_gene: test_gene.o 
	$(CC) $(CFLAGS) -c $(V0)/tests/test_gene.c -o ./output/test_gene.o $(LDFLAGS)	
	
gene: 
	$(CC) $(CFLAGS) -c $(V0)/src/gene.c -o ./output/gene.o $(LDFLAGS)	

run_test_gene: test_gene
	./$(BIN)test_gene &

DNA : 
	python3 $(V0)/python/setup.py build
	cp build/lib*/*.so $(BIN)

test_DNA : 
	python3 -m pytest -s $(V0)/tests/test_DNA.py


# Binary optimized library
test_gene_bin.o: gene_bin


test_gene_bin: test_gene_bin.o 
	$(CC) $(CFLAGS) -c $(V1)/tests/test_gene_bin.c -o ./output/test_gene_bin.o $(LDFLAGS)	
	
gene_bin: 
	$(CC) $(CFLAGS) -c $(V1)/src/gene_bin.c -o ./output/gene_bin.o $(LDFLAGS)	

run_test_gene_bin: test_gene_bin
	./$(BIN)test_gene_bin &

DNA_bin : 
	python3 $(V1)python/setup_bin.py build
	cp build/lib*/*.so $(BIN)

test_DNA_bin : 
	python3 -m pytest -s $(V1)/tests/test_DNA_bin.py
