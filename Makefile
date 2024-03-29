CC = gcc

CFLAGS = -g -std=c99 -Wall -Wextra -fopenmp
CTESTSFLAGS=-Wno-unused-function -Wno-unused-parameter
LDFLAGS = -lcmocka
MPIFLAGS=--mca opal_warn_on_missing_libcuda 0

V0 = ./versions/naive
V1 = ./versions/bin
V2 = ./versions/bool
V3 = ./versions/parallel_bin

# mpich or openmpi
MPI=.openmpi
NP=4

V0SRC=$(V0)/src
V1SRC=$(V1)/src
V2SRC=$(V2)/src
V3SRC=$(V3)/src


V0TESTS=$(V0)/tests
V1TESTS=$(V1)/tests
V2TESTS=$(V2)/tests
V3TESTS=$(V3)/tests


V0PYTHON=$(V0)/python
V1PYTHON=$(V1)/python
V2PYTHON=$(V2)/python
V3PYTHON=$(V3)/python


BIN = ./bin
OUTPUT = ./output
BUILD = ./build

$(shell mkdir -p $(BIN) $(OUTPUT) $(BUILD))

.PHONY: clean all check build install


ifeq (run, $(firstword $(MAKECMDGOALS)))
  runargs := $(wordlist 2, $(words $(MAKECMDGOALS)), $(MAKECMDGOALS))
  $(eval $(runargs):;@true)
endif

#For all compiling, building and executing
all: DNA check 

#For installing python setup and setup_bin
install: setup_naive setup_bin setup_bool

setup_naive:
	python3 $(V0PYTHON)/setup.py install

setup_bin:
	python3 $(V1PYTHON)/setup_bin.py install

setup_bool:
	python3 $(V2PYTHON)/setup_bool.py install	

#For only building and testing interface
build: DNA DNA_bin DNA_bool

#Clean compilation files & output result
clean :
	find . -type d  -name "__pycache__" -exec rm -rv {} +
	find . -type f  -name "*.so" -exec rm -rv {} +
	rm -rf $(BUILD) $(BIN)
	rm -rf $(V3)/build
	rm -rf $(V3PYTHON)/DNA_mod.c

#For only executing tests
check: run_test_gene test_DNA run_test_gene_bin test_DNA_bin test_DNA_bool

#For only running the non-binary program
run:
	python3 $(V0PYTHON)/setup.py install
	python3 $(V0PYTHON)/main.py $(runargs)


#For only running the binary program
run_bin:
	python3 $(V1PYTHON)/setup_bin.py install
	python3 $(V1PYTHON)/main_bin.py $(runargs)

setup_bin_par:
	sed "s|sources.*|sources=[\"$(V3PYTHON)/DNA_mod.pyx\"],|g" -i $(V3PYTHON)/setup_bin.py
	MPICC=mpicc python3 $(V3PYTHON)/setup_bin.py build_ext --inplace
	mv *.so $(V3PYTHON)

run_bin_par: setup_bin_par
	mpirun $(MPIFLAGS) -np $(NP) python3 $(V3PYTHON)/main.py $(runargs)

%.o:
	$(CC) $(CFLAGS) -c -o $(BUILD)/$(*F).o $(*D)/$(*F).c


# Naive library
.PHONY: gene

gene: $(V0SRC)/gene.o

test_gene: gene $(V0TESTS)/test_gene.o
	$(CC) $(CFLAGS) -o $(BIN)/$@ $(BUILD)/$@.o $(LDFLAGS)	

run_test_gene: test_gene
	$(BIN)/test_gene &

DNA : 
	python3 $(V0PYTHON)/setup.py build
	cp build/lib*/*.so $(BIN)

test_DNA : 
	python3 -m pytest -s $(V0TESTS)/test_DNA.py


# Binary optimized library
.PHONY: gene_bin

gene_bin: $(V1SRC)/gene_bin.o

test_gene_bin: gene_bin $(V1TESTS)/test_gene_bin.o
	$(CC) $(CFLAGS) $(CTESTSFLAGS) -o $(BIN)/$@ $(BUILD)/$@.o $(LDFLAGS)

run_test_gene_bin: test_gene_bin
	$(BIN)/test_gene_bin &

DNA_bin : 
	python3 $(V1PYTHON)/setup_bin.py build
	cp build/lib*/*.so $(BIN)

test_DNA_bin : 
	python3 -m pytest -s $(V1TESTS)/test_DNA_bin.py
	
# Boolean version

.PHONY: gene_bool

gene_bool: $(V2SRC)/gene_bool.o

#test_gene_bool: gene_bool $(V2TESTS)/test_gene_bool.o
#	$(CC) $(CFLAGS) -o $(BIN)/$@ $(BUILD)/$@.o $(LDFLAGS)

run_test_gene_bool: test_gene_bool
	$(BIN)/test_gene_bool &

DNA_bool : 
	python3 $(V2PYTHON)/setup_bool.py build
	cp build/lib*/*.so $(BIN)

test_DNA_bool : 
	python3 -m pytest -s $(V2TESTS)/test_DNA_bool.py
