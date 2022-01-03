from distutils.core import setup, Extension

DNAb_module = Extension("DNA_bin", sources = [ "gene_bin.c", "DNA_bin.c" ])

setup(name        = "DNA_bin",
      version     = "2.0",
      description = "Genome analyser", 
      ext_modules = [ DNAb_module ])
