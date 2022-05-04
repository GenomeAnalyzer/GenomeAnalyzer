#!/bin/bash

gcc -c -o gene_bin versions/bin/src/gene_bin.c

gcc -march=native -Ofast -funroll-loops -finline-functions -ftree-vectorize -o intr test_intrinsic.c gene_bin

./intr

#rm intr gene_bin

#maqao oneview -R1 --compiler="GNU" --language="c" --uarch="KNM" --envv_OMP_NUM_THREADS="1" -- ./intr
