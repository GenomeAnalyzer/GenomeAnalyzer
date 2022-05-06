# About

This project is a genome analyzer, it allows you to obtain mRNA, mutations, acidic amino chains, and gene detection matching of .fasta sequences.

# Requirement 

- Python 3

- Cython ( it's a python module)

- MPI ( for the parallel version)

# Installation

`make setup_VERSION`

# Run

In order to run you need to have a .fasta file in the subfolder ./fastas/

All sub command to run the program:

|  Long command   | Short command |                 function                  |
| :-------------: | :-----------: | :---------------------------------------: |
|     --help      |      -h       |         display those information         |
|    --output     |      -o       |     output information in html files      |
| --sequences <n> |    -s <n>     | sequences to be analyzed (without = max ) |

`make run_VERSION <sub command>` 

# Output

This project output the result in a HTML format into the ./output/ folder.

# Versions

There are 3 main versions:

- Bin: it's a binary version (an array of int that contains 32 bits)

- Naive: A naive version with only 1 bit per cell on the array

- Parallel bin: It's the binary version but with parallelism( MPI & OPENMP) and intrinsic 