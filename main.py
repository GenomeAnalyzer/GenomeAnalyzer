#coding: utf8
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
    #print(str(sequence[i]))

    gene.append(DNA.detecting_genes(array.array('I',sequence[i])))
    print(str(len(gene[i])))



    for j in range (len(gene[i])-1):
        res = (DNA.generating_mRNA(array.array('i',sequence[i][gene[i][j][0]:(gene[i][j][1])])))
        if res:
            mrna.append(res.decode('cp1252'))
        print("sequence = "+str(sequence[i][gene[i][j][0]:(gene[i][j][1])]))
        res =  DNA.generating_amino_acid_chain(array.array('i',sequence[i][gene[i][j][0]:(gene[i][j][1])]))
        print('résultat python ='+str(res.decode('cp1252')))

        if res:

            chain.append(res.decode('cp1252'))



    f.write("\nMRNA:")
    f.write(str(mrna))
    f.write("\nCHAIN:")
    f.write(str(chain))    

    f.write("\nMatching:")
    if (len(gene[i])-1 > 2):
        for j in range((int((len(gene[i])-1)/2)), -1, -1):
            for k in range((int((len(gene[i])-1)/2))):
                    f.write("["+str(j)+","+str(k)+"]=")
                    matching.append((DNA.calculating_matching_score(array.array('H',gene[i][j]),array.array('H',gene[i][k]))))
                    f.write(str(matching[k])+"\t")
    i=i+1
f.close()