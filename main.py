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

matching = [[]]
#gene = [[]]

f = open("rapport.txt", "a")

for file in glob.glob("fastas/*.fasta"):
    f.write("######################"+file+"######################")

    read = m.read_file(file)

    sequence.append(DNA.convert_to_binary(read,len(read)))

    gene.append(DNA.detecting_genes(array.array('I',sequence[i])))

    for j in range (len(gene[i])-1):
         mrna.append((DNA.generating_mRNA(array.array('i',sequence[i][gene[i][j][0]:(gene[i][j][1])]))))
         #chain[i].append(DNA.generating_amino_acid_chain(array.array('i',sequence[i][gene[i][j][0]:(gene[i][j][1])])))
    f.write("\nMRNA:")
    f.write(str(mrna))

    f.write("\nMatching:")
    for j in range((int((len(gene[i])-1)/2)), -1, -1):
        for k in range((int((len(gene[i])-1)/2))):
                f.write("["+str(j)+","+str(k)+"]=")
                matching.append((DNA.calculating_matching_score(array.array('H',gene[i][j]),array.array('H',gene[i][k]))))
                f.write(str(matching[k])+"\t")
    i=i+1
f.close()