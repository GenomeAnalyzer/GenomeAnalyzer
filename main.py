import DNA
import array
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

    read = m.read_file(file)

    sequence[i] = DNA.convert_to_binary(read,len(read))

    gene[i] = DNA.detecting_genes(array.array('I',sequence[i]))

    for j in range (len(gene[i])-1):
         mrna[i].append(DNA.generating_mRNA(array.array('i',sequence[i][gene[i][j][0]:(gene[i][j][1])])))
         chain[i].append(DNA.generating_amino_acid_chain(array.array('i',sequence[i][gene[i][j][0]:(gene[i][j][1])])))
       


    for j in range((int((len(gene[i])-1)/2)), -1, -1):
        for k in range((int((len(gene[i])-1)/2))):
                matching[i].append(DNA.calculating_matching_score(array.array('H',gene[i][j]),array.array('H',gene[i][k])))


    i+=1
