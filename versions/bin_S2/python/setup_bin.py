### From a setup.py on a github repo, we don't own the code

import os
import sys
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

extra_flags_c = ["-fopenmp","-O3", "-march=native" ]
extra_link_args_c = ["-fopenmp" ,"-O3", "-march=native"]
mpi_compile_args = []
mpi_link_args = []

# If true, do not check for parallel
avoid_parallel_test = False

# Get the MPI from environmental variables
parallel = False
if "MPICC" in os.environ:
    mpicc = os.environ["MPICC"]
    parallel = True
    print()
    print("Parallel compiler setted to:", mpicc)
    print()

# Check for the python parallel libraries
python_parallel = True
try:
    import pypar
except:
    try:
        import mpi4py
    except:
        #parallel = False
        python_parallel = False

# Setup the parallel environemnt
if parallel:
    # If we are here we can compile using MPI support
    mpi_compile_args = os.popen("%s -show" %
                                mpicc).read().strip().split(' ')[1:]
    mpi_link_args = os.popen("%s -show" % mpicc).read().strip().split(' ')[1:]
    extra_flags_c += ["-D_MPI"]

print(mpi_compile_args)
# Check if it is python2 or 3
if sys.version_info[0] < 3:
    print("Running on python 2, added the flag -D_PYTHON2")
    extra_flags_c += ["-D_PYTHON2"]


ext_modules = [
    Extension("DNA_mod",
              sources=["./versions/bin_S2/python/DNA_mod.pyx"],
              language='c',
              extra_compile_args=extra_flags_c + mpi_compile_args,
              extra_link_args=mpi_link_args + extra_link_args_c
              )
]

setup(
    name="DNA_mod",
    cmdclass={"build_ext": build_ext},
    ext_modules=ext_modules
)

if not python_parallel and not parallel and not avoid_parallel_test:
    print()
    print("======= WARNING =======")
    print("Nor python parallel neither MPI compiler found.")
    print("If you whish to activate MPI acceleration,")
    print("Consider installing either pypar or mpi4py")
    print("For example, try to run: ")
    print(" >>> MPICC=mpicc python3 " + " ".join(sys.argv))
    print("Note: clean the build directory if you whish to recompile the code.")
    print("=======================")
    print()
elif not parallel and not avoid_parallel_test:
    print()
    print("======= WARNING =======")
    print("No MPI compiler found, please specify MPICC environmental variable")
    print("For example, try to run: ")
    print(" >>> MPICC=mpicc python3 " + " ".join(sys.argv))
    print("Note: clean the build directory if you whish to recompile the code.")
    print("=======================")
    print()
elif not python_parallel and not avoid_parallel_test:
    print()
    print("======= WARNING =======")
    print("No Python MPI library found")
    print("Supported libraries:")
    print(" - pypar ")
    print(" - mpi4py ")
    print()
    print("Note: Fast MPI implemetation will crash if used")
    print("      consider to install one of these libraries.")
    print("      (No need to reinstall python-sscha)")
    print("=======================")
    print()
elif not avoid_parallel_test:
    print()
    print(" PARALLEL ENVIRONMENT DETECTED CORRECTLY ")
    print()
