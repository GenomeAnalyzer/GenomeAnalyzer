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
all: gcc_bit  clang_bit  icc_bit gcc_st clang_st icc_st
else
all:  gcc_bit  llvm_bit gcc_st llvm_st
endif


gcc_plot: main_plot.c ../gene_st.c
	$(GCC) $(GCCFLAGS) -g -o $@ $^

gcc_bit: main_bit.c ../gene_bit.c
	$(GCC) $(GCCFLAGS) -g -o $@ $^

llvm_bit: main_bit.c ../gene_bit.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^

clang_bit: main_bit.c ../gene_bit.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^

icc_bit: main_bit.c ../gene_bit.c
	$(ICC) $(ICCFLAGS) -o $@ $^

gcc_st: main_st.c ../gene_st.c
	$(GCC) $(GCCFLAGS) -g -o $@ $^

llvm_st: main_st.c ../gene_st.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^

clang_st: main_st.c ../gene_st.c
	$(CLANG) $(CLANGFLAGS) -o $@ $^

icc_st: main_st.c ../gene_st.c
	$(ICC) $(ICCFLAGS) -o $@ $^	

clean :
	@rm -f *.o  gcc_bit clang_bit  llvm_bit icc_bit gcc_st clang_st llvm_st icc_st
