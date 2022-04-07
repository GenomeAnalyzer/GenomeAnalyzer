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
DNAb_module = Extension("DNA_bin", sources = [ "./versions/bin/src/gene_bin.c", "./versions/bin/src/DNA_bin.c" ],extra_compile_args = extra_compile_args)

setup(name        = "DNA_bin",
      version     = "2.0",
      description = "Genome analyser", 
      ext_modules = [ DNAb_module ])
