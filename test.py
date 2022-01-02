#coding: utf8
import DNA_bin
import DNA
import array
import glob
import moduleDNA as m
import os


read = []
sequence = [[]]
sequence2 = [[]]

mrna = [[[]]]
i = 0
gene1 = list((()))

gene2 = list((()))


chain = [[[]]]

mutation = [[[]]]

matching = [[]]
#gene = [[]]

'''
sequence =DNA_bin.convert_to_binary("ATGAAAAAAAAAAAAAAAGTAA",len("ATGAGGTGTTAAAAAAAAAGTAA")) 

print(str(sequence))

#sequence[0] = int(bin(sequence[0])[:1:-1], 2)

print(str(bin(sequence[0])))

gene1 =DNA_bin.detecting_genes(array.array('l',sequence))
print('non')
print(str(gene1))

#print( " seq = "+str(DNA_bin.get_piece_binary_array(array.array('l',sequence[i]),gene1[i][0][0],gene1[i][0][1]- gene1[i][0][0])))


'''
for file in glob.glob("fastas/*.fasta"):

    read = m.read_file(file)

    sequence.append(DNA_bin.convert_to_binary(read,len(read)))

    sequence2.append(DNA.convert_to_binary(read,len(read)))

    gene1.append(DNA_bin.detecting_genes(array.array('l',sequence[i])))

    gen1 = DNA_bin.detecting_genes(array.array('l',sequence[i]))

    gene2.append(DNA.detecting_genes(array.array('L',sequence2[i])))

    gen2 = DNA.detecting_genes(array.array('L',sequence2[i]))

    #print("gene2"+str(gene2))
    print(str(file))

    for j in range (len(gene2[i])-1):
        print('####')
        #print("gene2"+str(gene1))

        print("taille = "+str(len(gen2))+" "+str(len(gen1)))
        print(""+str(gene2[i][j][0]) +" "+str(gene2[i][j][1]))
        print(""+str(gene1[i][j][0]) +" "+str(gene1[i][j][1]))

        print(" res = "+str(sequence2[i][gene2[i][j][0]:gene2[i][j][1]]))

        #print("seq = "+str(m.get_sequence(sequence[i],gene1[i][j][0],(gene1[i][j][1]))))

        print( " seq = "+str(bin(DNA_bin.get_piece_binary_array(array.array('l',sequence[i]),gene1[i][j][0],gene1[i][j][1]- gene1[i][j][0])[0])))

        if (len(gene1[i])-1 > 2):
            for j in range((int((len(gene1[i]))/2)), -1, -1):
                for k in range((int((len(gene1[i]))/2))):     
                    res = (DNA.calculating_matching_score(array.array('H',sequence2[i][gene2[i][j][0]:gene2[i][j][1]]),array.array('H',sequence2[i][gene2[i][k][0]:gene2[i][k][1]])))
                    print(" matching = "+str(res))
                    res2 = DNA_bin.calculating_matching_score(array.array('l',sequence[i]),gene1[i][j][0],gene1[i][j][1]- gene1[i][j][0],array.array('l',sequence[i]),gene1[i][k][0],gene1[i][k][1]- gene1[i][k][0])

                    print(" matching 2 ="+str(res2))

                    print('#####')


    i = i+1

    if i ==4:
        break


'''
        print(" mrna = "+str(DNA_bin.generating_mRNA(array.array('l',sequence[i]),gene1[i][j][0],gene1[i][j][1]- gene1[i][j][0])))

        print(" mrna 2"+ str(DNA.generating_mRNA(array.array('H',sequence2[i][gene2[i][j][0]:gene2[i][j][1]+1]))))
        print('test')

        print(" amino = "+str(DNA_bin.generating_amino_acid_chain(array.array('l',sequence[i]),gene1[i][j][0],gene1[i][j][1]- gene1[i][j][0])))

        print(" amino 2="+ str(DNA.generating_amino_acid_chain(array.array('H',sequence2[i][gene2[i][j][0]:gene2[i][j][1]+1]))))

        print(" mutation = "+str(DNA_bin.detecting_mutations(array.array('l',sequence[i]),gene1[i][j][0],gene1[i][j][1]- gene1[i][j][0])))

        print(" mutation 2="+ str(DNA.detecting_mutations(array.array('H',sequence2[i][gene2[i][j][0]:gene2[i][j][1]+1]))))        
'''