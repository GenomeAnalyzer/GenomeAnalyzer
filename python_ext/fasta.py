import FASTAOP
import os
import json

fastas_folder = "../fastas"
seq_filename = f"{fastas_folder}/NC_045512.2.fasta"

if not os.path.exists(seq_filename):
    print("File does not exists")
    exit()


# Load codons
fo = open("../codons.json","r")
codons = json.load(fo)
fo.close()
CODON_START = "AUG"
CODON_STOPS = ["UAA","UAG","UGA"]


fo = open(seq_filename)
lines = fo.readlines()

seq = ""

for line in lines[1:]:
    line = line[:-1]
    seq += line

seq2 = seq.replace("T","U")

print(seq2)

seqs = [[],[],[]]

for i in range(len(seq2) - 2):
    c = seq2[i]
    codon = seq2[i:i+3]
    if codon == CODON_START:
        res = f"{i} Sta"
        j = i + 3

        n_seq = i%3

        while j < len(seq2) - 2:
            codon = seq2[j:j+3]
            if not codon in codons: print("invalid codon"); input()
            else: res += f" {codons[codon]}"
            if codon in CODON_STOPS:
                print(res)
                print((len(res)+1)/4)
                # input()
                break
            j += 3

        seqs[n_seq].append(res)

for i,s in enumerate(seqs):
    # print(len(s))
    n = 0
    for ss in s:
        n += len(ss)
    print(i, len(s), n)