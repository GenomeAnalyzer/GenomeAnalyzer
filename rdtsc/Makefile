IFICC = $(shell icc --version | head -c3)
$(info if_icc=$(IFICC))

MACHINE = $(shell uname -n)
$(info machine=$(MACHINE))

TYPECLANG = $(shell clang --version | head -c16)
$(info type_clang=$(TYPECLANG))

GCC = gcc
CLANG = clang
ICC = icc

ifeq ($(MACHINE),snb01)
	GCCFLAGS = -g -march=sandybridge -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize
	CLANGFLAGS = -g -march=sandybridge -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize
	ICCFLAGS = -g -xhost -mavx2 -Ofast -funroll-all-loops -finline-functions
else 
	GCCFLAGS = -g -march=native -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize
	CLANGFLAGS = -g -march=native -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize
	ICCFLAGS = -g -xhost -mavx2 -Ofast -funroll-all-loops -finline-functions
endif 

.PHONY: clean all check

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

ifeq ($(IFICC),icc)
all: gcc_nobin gcc_bin gcc_bool clang_nobin clang_bin clang_bool icc_nobin icc_bin icc_bool
else
all: gcc_nobin gcc_bin gcc_bool llvm_nobin llvm_bin llvm_bool
endif

gcc_nobin: main.c ../versions/naive/src/gene.c
	$(GCC) $(GCCFLAGS) -o $@ $^
	
gcc_bin: main_bin.c ../versions/bin/src/gene_bin.c
	$(GCC) $(GCCFLAGS) -o $@ $^

gcc_bool: main_bool.c ../versions/bool/src/gene_bool.c
	$(GCC) $(GCCFLAGS) -o $@ $^

llvm_nobin: main.c ../versions/naive/src/gene.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^
	
llvm_bin: main_bin.c ../versions/bin/src/gene_bin.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^

llvm_bool: main_bool.c ../versions/bool/src/gene_bool.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^

clang_nobin: main.c ../versions/naive/src/gene.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^
	
clang_bin: main_bin.c ../versions/bin/src/gene_bin.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^

clang_bool: main_bin.c ../versions/bool/src/gene_bool.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^

icc_nobin: main.c ../versions/naive/src/gene.c
	$(ICC) $(ICCFLAGS) -o $@ $^
	
icc_bin: main_bin.c ../versions/bin/src/gene_bin.c
	$(ICC) $(ICCFLAGS) -o $@ $^

icc_bool: main_bin.c ../versions/bool/src/gene_bool.c
	$(ICC) $(ICCFLAGS) -o $@ $^

clean :
	@rm -f *.o gcc_nobin gcc_bin gcc_bool clang_nobin clang_bin clang_bool llvm_nobin llvm_bin llvm_bool icc_nobin icc_bin icc_bool
