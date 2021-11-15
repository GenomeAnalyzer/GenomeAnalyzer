from distutils.core import setup, Extension

DNA_module = Extension("DNA", sources = [ "gene.c", "DNA.c" ])

setup(name        = "DNA",
      version     = "0.1",
      description = "Genome analyser", 
      ext_modules = [ DNA_module ])