
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
mrna = [[[]]]
i = 0
gene = [[[]]]

chain = [[[]]]

mutation = [[[]]]

matching = [[[]]]
#gene = [[]]
for file in glob.glob("fastas/*.fasta"):

    #print(str(file))
    read = m.read_file(file)
    #print(test)
#    print(str(len(test)))
    sequence[i] = DNA.convert_to_binary(read,len(read))
    #print(str(sequence))

    gene[i] = DNA.detecting_genes(array.array('I',sequence[i]))
    #print (str(gene[i][0][0])+" "+str(gene[i][0][1])+ " " + str(len(sequence[i])))

    for j in range (len(gene[i])-1):
         mrna[i].append(DNA.generating_mRNA(array.array('i',sequence[i][gene[i][j][0]:(gene[i][j][1])])))
         chain[i].append(DNA.generating_amino_acid_chain(array.array('i',sequence[i][gene[i][j][0]:(gene[i][j][1])])))
    print(chain)
       


    for j in range((int((len(gene[i])-1)/2)), -1, -1):
        for k in range((int((len(gene[i])-1)/2))):
                matching[i].append(DNA.calculating_matching_score(array.array('H',gene[i][j]),array.array('H',gene[i][k])))
    #print('matching = '+str(matching))




        #mrna[i][j]= DNA.generating_mRNA(array.array('i',sequence[i][gene[i][j][0]:(gene[i][j][1])]))
        #bytes(s, 'utf8')

       # mutation[i][j] = DNA.detecting_mutations(


    #print()

    i+=1
    break 
    #ok = np.ctypeslib.as_array(truc, shape=(3072))
#print((str(gene)))