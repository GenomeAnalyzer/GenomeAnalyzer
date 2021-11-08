from distutils.core import setup, Extension

FASTAOP_module = Extension("FASTAOP", sources = [ "FASTAOP_core.c", "FASTAOP.c" ])

setup(name        = "FASTAOP",
      version     = "0.1",
      description = "FASTA OPerations library", 
      license     = "MIT",
      ext_modules = [ FASTAOP_module ])