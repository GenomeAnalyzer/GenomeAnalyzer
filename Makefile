CC = gcc

CFLAGS = -g -std=c99 -Wall

LDFLAGS = -lcmocka

V0 = ./versions/naive
V1 = ./versions/bin
V2 = ./versions/bool

V0SRC=$(V0)/src
V1SRC=$(V1)/src
V2SRC=$(V2)/src

V0TESTS=$(V0)/tests
V1TESTS=$(V1)/tests
V2TESTS=$(V2)/tests


V0PYTHON=$(V0)/python
V1PYTHON=$(V1)/python
V2PYTHON=$(V2)/python

BIN = ./bin
OUTPUT = ./output
BUILD = ./build

$(shell mkdir -p $(BIN) $(OUTPUT))

.PHONY: clean all check


ifeq (run, $(firstword $(MAKECMDGOALS)))
  runargs := $(wordlist 2, $(words $(MAKECMDGOALS)), $(MAKECMDGOALS))
  $(eval $(runargs):;@true)
endif

#For all compiling, building and executing
all: DNA check 

#For installing python setup and setup_bin
install:
	sudo python3 $(V0PYTHON)/setup.py install
	sudo python3 $(V1PYTHON)/setup_bin.py install
	sudo python3 $(V2PYTHON)/setup_bool.py install	

#For only building and testing interface
build: DNA DNA_bin DNA_bool

#Clean compilation files & output result
clean :
	find . -type d  -name "__pycache__" -exec rm -rv {} +
	rm -rf  build .pytest_cache *.so
	rm -rf $(BUILD) $(BIN)/*

#For only executing tests
check: run_test_gene test_DNA run_test_gene_bin test_DNA_bin test_DNA_bool

#For only running the non-binary program
run:
	sudo python3 $(V0PYTHON)/setup.py install
	python3 $(V0PYTHON)/main.py $(runargs)
	

#For only running the binary program
run_bin:
	sudo python3 $(V1PYTHON)/setup_bin.py install
	python3 $(V1PYTHON)/main_bin.py $(runargs)


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
	$(CC) $(CFLAGS) -o $(BIN)/$@ $(BUILD)/$@.o $(LDFLAGS)

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
