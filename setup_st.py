from distutils.core import setup, Extension

DNAst_module = Extension("DNA_st", sources = [ "gene_st.c", "DNA_st.c" ])

setup(name        = "DNA_bit",
      version     = "2.0",
      description = "Genome analyser", 
      ext_modules = [ DNAst_module ])
