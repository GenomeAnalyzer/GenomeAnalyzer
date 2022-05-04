#cython: language_level=3

cdef extern from "gene_bin.c":
   void launch()

def py_launch():
    launch()