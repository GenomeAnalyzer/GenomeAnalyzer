CC = gcc

CFLAGS = -g -std=c99 -Wall -Wextra -fopenmp
CTESTSFLAGS=-Wno-unused-function -Wno-unused-parameter
LDFLAGS = -lcmocka
MPIFLAGS=--mca opal_warn_on_missing_libcuda 0

SRC=src
TESTS=tests
PYTHON=python

BIN = .
OUTPUT = ./output
BUILD = ./build

# mpich or openmpi
MPI=.openmpi
NP=4


$(shell mkdir -p $(OUTPUT) $(BUILD))

.PHONY: clean all check build install


ifeq (run, $(firstword $(MAKECMDGOALS)))
  runargs := $(wordlist 2, $(words $(MAKECMDGOALS)), $(MAKECMDGOALS))
  $(eval $(runargs):;@true)
endif

#For all compiling, building and executing
all: install check

#For installing python setup and setup_bin
install: setup_bin_par

setup_bin_par:
	sed "s|sources.*|sources=[\"$(PYTHON)/DNA_mod.pyx\"],|g" -i $(PYTHON)/setup_bin.py
	MPICC=mpicc python3 $(PYTHON)/setup_bin.py build_ext --inplace

#Run the parallel binary version
run_bin_par: setup_bin_par
	mv *.so  $(PYTHON)
	mpirun $(MPIFLAGS) -np $(NP) python3 $(PYTHON)/main.py $(runargs)

# /usr/lib/python3.10/site-packages/DNA_mod.cpython-310-x86_64-linux-gnu.so

%.o:
	$(CC) $(CFLAGS) -c -o $(BUILD)/$(*F).o $(*D)/$(*F).c

gene_bin: $(SRC)/gene_bin.o

test_gene_bin: gene_bin $(TESTS)/test_gene_bin.o
	$(CC) $(CFLAGS) $(CTESTSFLAGS) -o $(BIN)/$@ $(BUILD)/$@.o $(LDFLAGS)

run_test_gene_bin: test_gene_bin
	$(BIN)/test_gene_bin

test_DNA_bin : 
	python3 -m pytest -s $(TESTS)/test_DNA_bin.py

clean:
	rm -rf $(BUILD) $(PYTHON)/*.so
	rm -rf $(PYTHON)/DNA_mod.c