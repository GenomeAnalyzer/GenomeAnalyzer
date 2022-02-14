import numpy as np
import DNA_bin
import moduleDNA as m
import array
import os
from math import *
import time

def convertbin(seq_char):
	seq_bin_char = seq_char.replace("A", "00")
	seq_bin_char = seq_bin_char.replace("G", "01")
	seq_bin_char = seq_bin_char.replace("C", "10")
	seq_bin_char = seq_bin_char.replace("T", "11")
	seq_bin_list = seq_bin_char.split()

	seq_bin_list = [np.int64(np.uint64(int(i,2))) for i in seq_bin_list]
	
	return seq_bin_list

file = "fastas/LC528232.1.fasta.txt"
read = m.read_fasta_txt(file)

# print(read)
size = ceil(len(read.replace('\n',''))/32.0)

start = time.monotonic()
seq_bin = DNA_bin.convert_to_binary(read,len(read),size)
end = time.monotonic()
print("Converter time in C : ", end-start)

start = time.monotonic()
seq_bin2 = convertbin(read)
end = time.monotonic()
print("Converter time in Python : ", end-start)

for i in range(len(seq_bin2)):
	if(seq_bin[i] ^ seq_bin2[i]):
		print("NON EGAUX")

readn = read.replace('\n','')
genes = DNA_bin.detecting_genes(array.array('l',seq_bin))
for i in genes:
	i[0] = int(i[0] >> 1)
	i[1] = int(i[1] >> 1)
	print(readn[i[0]:i[1]+1])
