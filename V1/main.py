from os import read
import DNA
import array
import numpy as np
import faulthandler
import sys
import ctypes 
import glob
import moduleDNA as m


read = []
sequence = [[]]
mrna = [[]]
i = 0
gene = [[[]]]
#gene = [[]]
for file in glob.glob("fastas/*.fasta"):

    #print(str(file))
    read = m.read_file(file)
    #print(test)
#    print(str(len(test)))
    sequence[i] = DNA.convert_to_binary(read,len(read))
    #print(str(sequence))

    gene[i] = DNA.detecting_genes(array.array('I',sequence[i]))
    print (str(gene[i][0][0])+" "+str(gene[i][0][1])+ " " + str(len(sequence[i])))
    for j in range(len(gene[i])):
        mrna.append(DNA.generating_mRNA(array.array('i',sequence[i][gene[i][j][0]:(gene[i][j][1])])))

        print(mrna, end=' ')
        print()

    i+=1
    break 
    #ok = np.ctypeslib.as_array(truc, shape=(3072))
print((str(gene)))