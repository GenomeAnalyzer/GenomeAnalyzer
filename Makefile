CC = gcc

CFLAGS = -g -std=c99

LDFLAGS = -lcmocka

.PHONY: clean all check

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

#For all compiling, building and executing
all: DNA check 

#For only building and testing interface
build: DNA DNA_bin

#For only executing tests
check: test_gene_bin test_gene test_DNA test_DNA_bin

#For only running the non-binary program
run:
	sudo python3 setup.py install
	python3 main.py

#For only running the binary program
run_bin:
	sudo python3 setup_bin.py install
	python3 main_bin.py

test_gene: test_gene.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_gene.o: gene.c

test_gene_bin: test_gene_bin.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_gene_bin.o: gene_bin.c

DNA : 
	python3 setup.py build
	cp build/lib*/*.so .

DNA_bin : 
	python3 setup_bin.py build
	cp build/lib*/*.so .

test_DNA : 
	python3 -m pytest -s

test_DNA_bin : 
	python3 -m pytest -s

clean :
	@rm -rf __pycache__ build .pytest_cache output *.so
	@rm -f *.o test_gene test_gene_bin

