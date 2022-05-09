#cython: language_level=3

cdef extern from "../src/gene_bin.c":
   void launch(int output,int sequence)

def py_launch(int output,int sequence):
    launch(output,sequence)