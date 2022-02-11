import glob
#import os
from math import ceil

def read_file(name):
    gene = ""
    with open(name, "r") as file:
        rows = file.readlines()
        for row in rows[1:]:
            gene = gene + row.strip()
    return gene

def write_file(namefile, gene):
    l = ceil(len(gene)/32)
    with open(namefile, "w") as file:
        for i in range(0,l):
            file.write(gene[i*32:(i+1)*32]+'\n')

# if not os.path.exists('fasta_text'):
#     os.makedirs('fasta_text')

for file in glob.glob("fastas/*.fasta"):
    read = read_file(file)
    write_file(file+".txt", read)
    # write_file("fasta_text/"+file[7:-5]+"txt", read)
