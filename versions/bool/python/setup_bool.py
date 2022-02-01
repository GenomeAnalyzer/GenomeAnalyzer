from distutils.core import setup, Extension

DNAbi_module = Extension("DNA_bool", sources = [ "./versions/bool/src/gene_bool.c", "./versions/bool/src/DNA_bool.c" ])

setup(name        = "DNA_bool",
      version     = "2.0",
      description = "Genome analyser", 
      ext_modules = [ DNAbi_module ])
