from distutils.core import setup, Extension

DNA_module = Extension("DNA", sources = [ "./versions/naive/src/gene.c", "./versions/naive/src/DNA.c" ])

setup(name        = "DNA",
      version     = "1.0",
      description = "Genome analyser", 
      ext_modules = [ DNA_module ])
