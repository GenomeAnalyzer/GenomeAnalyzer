# About

This projet is a genome analyzer, it allows you to obtain mRNA, mutations, acidic amino chains, and gene detection matching of .fasta sequences.


# Installation

`make build`

# Run

In order to run you need to have a .fasta file in the subfolder ./fastas/

All sub command to run the program:

| Long command | Short command | function |
| :---: | :---: | :---: |
| --help | -h | display those information |
| --output | -o | output information in html files |
| --sequences <n> | -s <n> | sequences to be analyzed (without = max ) |

`make run <sub command>` 
# Output

This project output the result in an html format into the ./output/ folder.