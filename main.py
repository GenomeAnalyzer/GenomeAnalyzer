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

f = open("rapport.txt", "w")

for file in glob.glob("fastas/*.fasta"):
    f.write("\n######################"+file+"######################")

    read = m.read_file(file)

    sequence.append(DNA.convert_to_binary(read,len(read)))
    #print(str(sequence[i]))

    gene.append(DNA.detecting_genes(array.array('L',sequence[i])))
    #gene.append(DNA.detecting_genes(array.array('H',sequence[i])))

    f.write("\nMRNA\t\t\t\t\t\t\t\tCHAIN\t\t\t\t\t\t\t\tMUTATION\n")
    for j in range (len(gene[i])-1):
        #print('gene')       
        res = (DNA.generating_mRNA(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1])])))
        if res:
            mrna.append(res.decode("cp1252", "replace"))
            f.write(str(res.decode("cp1252", "replace")))
        #print('aminoacid ='+str(sequence[i][gene[i][j][0]:(gene[i][j][1])]))
        res2 =  DNA.generating_amino_acid_chain(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1])]))

        if res2:
            chain.append(res2.decode("cp1252", "replace"))
            f.write("\t\t\t\t\t\t\t\t"+str(res2.decode("cp1252", "replace")))
        #print('ok')
        #print('mutation =')
        mutres = DNA.detecting_mutations(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1])]))
        mutation.append(mutres)
        f.write("\t\t\t\t\t\t\t\t"+str(mutres)+"\n")




    #f.write("\nMRNA:")
    #f.write(str(mrna))
    #f.write("\nCHAIN:")
    #f.write(str(chain))    
    #f.write("\nMUTATION:")
    #f.write(str(mutation))

    f.write("\nMatching:")
    if (len(gene[i])-1 > 2):
        for j in range((int((len(gene[i])-1)/2)), -1, -1):
            for k in range((int((len(gene[i])-1)/2))):     
                        f.write("["+str(j)+","+str(k)+"]=")
                        matching.append((DNA.calculating_matching_score(array.array('H',gene[i][j]),array.array('H',gene[i][k]))))
                        f.write(str(matching[k])+"\t")
    i=i+1
f.close()