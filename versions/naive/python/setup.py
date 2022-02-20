from distutils.core import setup, Extension
import os
import sysconfig

_DEBUG = False
_DEBUG_LEVEL = 0

# Common flags for both release and debug builds.
extra_compile_args = sysconfig.get_config_var('CFLAGS').split()
extra_compile_args += [ "-Wall", "-Wextra"]
if _DEBUG:
    extra_compile_args += ["-g3", "-O0", "-DDEBUG=%s" % _DEBUG_LEVEL, "-UNDEBUG"]
else:
    extra_compile_args += ["-DNDEBUG", "-O3"]

#TODO add ,"-fopenmp"

DNA_module = Extension("DNA", sources = [ "./versions/naive/src/gene.c", "./versions/naive/src/DNA.c" ],extra_compile_args = extra_compile_args)

setup(name        = "DNA",
      version     = "1.0",
      description = "Genome analyser", 
      ext_modules = [ DNA_module ])
