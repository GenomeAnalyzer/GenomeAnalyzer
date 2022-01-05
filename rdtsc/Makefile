GCC = gcc
CLANG = clang
#ICC = icc

GCCFLAGS = -g -march=native -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize
CLANGFLAGS = -g -march=native -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize
#ICCFLAGS = -g -march=native -mavx2 -Ofast -vec -funroll-all-loops -finline-functions -xhost


.PHONY: clean all check

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

all: gcc_nobin gcc_bin clang_nobin clang_bin #icc_nobin icc_bin

gcc_nobin: main.c ../gene.c
	$(GCC) $(GCCFLAGS) -o $@ $^
	
gcc_bin: main_bin.c ../gene_bin.c
	$(GCC) $(GCCFLAGS) -o $@ $^

clang_nobin: main.c ../gene.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^
	
clang_bin: main_bin.c ../gene_bin.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^

icc_nobin: main.c ../gene.c
	$(ICC) $(ICCFLAGS) -o $@ $^
	
icc_bin: main_bin.c ../gene_bin.c
	$(ICC) $(ICCFLAGS) -o $@ $^

clean :
	@rm -f *.o gcc_nobin gcc_bin clang_nobin clang_bin #icc_nobin icc_bin