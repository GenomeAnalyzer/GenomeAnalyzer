#!/bin/bash

gcc -c -o gene_bin versions/bin/src/gene_bin.c

gcc -march=knl -mavx2 -Ofast -funroll-loops -finline-functions -ftree-vectorize -o intr test_intrinsic.c gene_bin

./intr

rm intr gene_bin
