# About

This projet is a genome analyzer, it allows you to obtain mRNA, mutations, acidic amino chains, and gene detection matching of .fasta sequences.

# Conversion ONLY for binary version

First of all, run this for converting .fasta files in .txt files with 32 letters per line: `make converter`

# Installation

For naive version : `make build DNA`
For binary version : `make build DNA_bin`
For boolean version : `make build DNA_bool`
For all versions : `make build`

# Tests

For the tests : `make check`

# Run

In order to run you need to have a .fasta file in the subfolder ./fastas/
x is the number of sequences to be analyzed
For naive version : `make run <x>`
For binary version : `make run_bin <x>`
For boolean version : `make run_bool <x>`

# Output

This project output the result in an html format into the ./output/ folder.