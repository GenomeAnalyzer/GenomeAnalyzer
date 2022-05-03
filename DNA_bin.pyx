#cython: language_level=3
cimport mpi4py.MPI as MPI
cimport mpi4py.libmpi as libmpi

cdef extern from "gene_bin.h":
   void launch(libmpi.MPI_Comm comm)

def py_launch(MPI.Comm comm):
    launch(comm.ob_mpi)