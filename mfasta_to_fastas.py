from os import mkdir
from os.path import exists
from sys import argv

# seq_filename = "sequences.fasta"
FASTA_CHAR = ">"
FASTAS_FILEPATH = "fastas"

# if not exists(seq_filename):
#     print(f"{seq_filename} : No such file")
#     exit()

if not exists(FASTAS_FILEPATH):
    mkdir(FASTAS_FILEPATH)

def mfasta_to_fastas(seq_filename):
    print(f"Unpacking {seq_filename} to {FASTAS_FILEPATH}/")
    fo = open(seq_filename)
    lines = fo.readlines()
    fo.close()
    if not lines[0].startswith(">"): 
        print(f"{seq_filename} : Incorrect file format, first line should start with the {FASTA_CHAR} char.\n")
        exit()

    for line in lines:
        if line.startswith(FASTA_CHAR):
            if "to_write" in locals():
                fw = open(filepath, "w")
                fw.write(to_write)
                fw.close()
            to_write = ""
            if not "|" in line : print("Incorrect line format :", line)
            splited = line.split("|")
            name = splited[0][1:].strip()
            # metadata = splited[1][:-1]
            filename = f"{name}.fasta"
            filepath = f"{FASTAS_FILEPATH}/{filename}"
        to_write += line
    fw = open(filepath, "w")
    fw.write(to_write)
    fw.close()

if __name__ == "__main__":
    if len(argv) < 2:
        print(f"No filepath specified\nUsage: python {argv[0]} [FASTA File]")
        exit()
    files = argv[1:]
    for file in files:
        if not file.endswith(".fasta"): print(f"{file} : Unknowned file format. Expected FASTA file format"); exit()
        if not exists(file): printf(f"{file} : No such file"); exit()
    
    for file in files: 
        mfasta_to_fastas(file)