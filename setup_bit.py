from distutils.core import setup, Extension

DNAbi_module = Extension("DNA_bit", sources = [ "gene_bit.c", "DNA_bit.c" ])

setup(name        = "DNA_bit",
      version     = "2.0",
      description = "Genome analyser", 
      ext_modules = [ DNAbi_module ])
