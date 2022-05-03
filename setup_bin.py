from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import os
import sysconfig

_DEBUG = False
_DEBUG_LEVEL = 0


mpi_compile_args = os.popen("mpicc --showme:compile").read().strip().split(' ')
mpi_link_args = os.popen("mpicc --showme:link").read().strip().split(' ')
# mpi_compile_args += ["-L /usr/lib/"]

#mpi_link_args = ["-L/usr/lib/ -pthread -lmpi -fopenmp"]

# Common flags for both release and debug builds.
extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
extra_compile_args += ["-Wall", "-Wextra"]
extra_compile_args += mpi_compile_args

if _DEBUG:
    extra_compile_args += ["-g3", "-O0", "-DDEBUG=%s" %
                           _DEBUG_LEVEL, "-UNDEBUG"]
else:
    extra_compile_args += ["-DNDEBUG", "-O3"]

here = os.path.abspath(os.path.dirname(__file__))
include_dirs = [here]
try:
    import mpi4py
except ImportError:
    MPI4PY = False
else:
    MPI4PY = True
    INCLUDE_MPI = '/usr/lib/openmpi/include'
    include_dirs.extend([
        INCLUDE_MPI,
        mpi4py.get_include()])

# TODO add ,"-fopenmp"
DNAb_module = Extension("DNA_bin", language='c', sources=["DNA_bin.pyx",
                        "gene_bin.c"],
                         extra_compile_args=extra_compile_args, extra_link_args=mpi_link_args)

setup(name="DNA_bin",
      version="2.0",
           cmdclass={"build_ext": build_ext},
      description="Genome analyser",
      ext_modules=[DNAb_module])
#cythonize(DNAb_module,compile_time_env={'MPI4PY': MPI4PY}))
if not MPI4PY:
    print('Warning: since importing mpi4py raises an ImportError,\n'
          '         the extensions are compiled without mpi and \n'
          '         will work only in sequencial.')
