
cimport mpi4py.MPI as MPI
cimport mpi4py.libmpi as libmpi

cdef extern from "./gene_bin.h":
    void launc() nogil
