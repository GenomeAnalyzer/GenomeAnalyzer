CC = gcc

CFLAGS = -g -std=c99 -Wall -Wextra -Wno-unused-function -Wno-unused-parameter
LDFLAGS = -lcmocka

SRC=src
TESTS=tests
PYTHON=python

BIN=.
OUTPUT = ./output
BUILD = ./build

$(shell mkdir -p $(BIN) $(OUTPUT) $(BUILD))

.PHONY: clean all check build install


ifeq (run, $(firstword $(MAKECMDGOALS)))
  runargs := $(wordlist 2, $(words $(MAKECMDGOALS)), $(MAKECMDGOALS))
  $(eval $(runargs):;@true)
endif

#For all compiling, building and executing
all: check 

#For installing python setup and setup_bin
install: setup_bin

setup_bin:
	sudo python3 $(PYTHON)/setup_bin.py install

#For only building and testing interface
build: DNA_bin

#Clean compilation files & output result
clean :
	find . -type d  -name "__pycache__" -exec rm -rv {} +
	rm -rf  build .pytest_cache *.so
	rm -rf $(BUILD)

#For only executing tests
check: run_test_gene_bin test_DNA_bin

#For only running the non-binary program
run: run_bin

#For only running the binary program
run_bin:
	sudo python3 $(PYTHON)/setup_bin.py install
	python3 $(PYTHON)/main_bin.py $(runargs)

%.o:
	$(CC) $(CFLAGS) -c -o $(BUILD)/$(*F).o $(*D)/$(*F).c

# Binary optimized library
.PHONY: gene_bin run_test_gene_bin DNA_bin test_DNA_bin

gene_bin: $(SRC)/gene_bin.o

test_gene_bin: gene_bin $(TESTS)/test_gene_bin.o
	$(CC) $(CFLAGS) -o $(BIN)/$@ $(BUILD)/$@.o $(LDFLAGS)

run_test_gene_bin: test_gene_bin
	$(BIN)/test_gene_bin &

DNA_bin : 
	python3 $(PYTHON)/setup_bin.py build
	cp build/lib*/*.so $(BIN)

test_DNA_bin : 
	python3 -m pytest -s $(TESTS)/test_DNA_bin.py
