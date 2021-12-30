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
fh = open('helloworld.html','w')
message = """<html>
<head><style> 
th, td {
  border-style: dotted;
  font-size: 10px; 
}
.title {
  border-style: dotted;
  font-size: 15px; 
} 
</style>
</head>

"""

fh.write(message)
for file in glob.glob("fastas/*.fasta"):
    message = "<details><summary>"+file+"</summary>"
    f.write("\n######################"+file+"######################")



    message += "<table>\n<tbody>\n<tr>\n<td class = \"title\">Sequence</td>\n<td class = \"title\">MRNA</td>\n<td class = \"title\">Chain</td>\n<td class = \"title\">Mutation</td>\n</tr>\n"
    

    read = m.read_file(file)

    sequence.append(DNA.convert_to_binary(read,len(read)))
    #print(str(sequence[i]))

    gene.append(DNA.detecting_genes(array.array('L',sequence[i])))
    #gene.append(DNA.detecting_genes(array.array('H',sequence[i])))

    f.write("\nMRNA\t\t\t\t\t\t\t\tCHAIN\t\t\t\t\t\t\t\tMUTATION\n")
    for j in range (len(gene[i])-1):
        message += "<tr><td>"+str(sequence[i][gene[i][j][0]:(gene[i][j][1])])+ "</td>\n"
        #print('gene')       
        res = (DNA.generating_mRNA(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1])])))
        if res:
            mrna.append(res.decode("cp1252", "replace"))
            f.write(str(res.decode("cp1252", "replace")))
            message+="<td>"+str(res.decode("cp1252", "replace"))+ "</td>\n"
        else:
            message+="<td> none </td>\n"

        #print('aminoacid ='+str(sequence[i][gene[i][j][0]:(gene[i][j][1])]))
        res2 =  DNA.generating_amino_acid_chain(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1])]))

        if res2:
            chain.append(res2.decode("cp1252", "replace"))
            f.write("\t\t\t\t\t\t\t\t"+str(res2.decode("cp1252", "replace")))
            message+="<td>"+str(res2.decode("cp1252", "replace"))+ "</td>\n"
        else:
            message+="<td> none </td>\n"            
        #print('ok')
        #print('mutation =')
        mutres = DNA.detecting_mutations(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1])]))
        mutation.append(mutres)
        message+="<td>"+str(mutres)+ "</td>\n</tr>\n"
        f.write("\t\t\t\t\t\t\t\t"+str(mutres)+"\n")
    message+= "</tbody>\n</table>\n"



    #f.write("\nMRNA:")
    #f.write(str(mrna))
    #f.write("\nCHAIN:")
    #f.write(str(chain))    
    #f.write("\nMUTATION:")
    #f.write(str(mutation))
    message += "<table>\n<tr>\n<th class = \"title\">Sequence</th>\n<th class = \"title\">Matching</th>\n</tr>\n<tbody>\n"

    f.write("\nMatching:")
    if (len(gene[i])-1 > 2):
        for j in range((int((len(gene[i]))/2)), -1, -1):
            for k in range((int((len(gene[i]))/2))):     
                message+="<tr><td>"+str(gene[i][j])+","+str(gene[i][k])+ "</td>\n"
                f.write("["+str(gene[i][j])+","+str(gene[i][k])+"]=")
                res = (DNA.calculating_matching_score(array.array('H',gene[i][j]),array.array('H',gene[i][k])))
                matching.append(res)
                message+="<td>"+str(res)+ "</td> \n</tr>\n"
                f.write(str(res)+"\t")
    message+= "</tbody>\n</table>\n</details>\n"  


    fh.write(message)
    i=i+1
message = "<table>\n<tr>\n<th class = \"title\">Sequence</th>\n<th class = \"title\">Matching</th>\n</tr>\n<tbody>\n"

f.write("\nMatching:")
for i in range((int((len(sequence))/2)), -1, -1):
    for c in range((int((len(sequence))/2))):     
        for j in (int((len(sequence[i])))):
            for k in (int((len(sequence[c])))):
                message+="<tr><td>"+str(gene[i][j])+","+str(gene[c][k])+ "</td>\n"
                f.write("["+str(gene[i][j])+","+str(gene[c][k])+"]=")
                res = (DNA.calculating_matching_score(array.array('H',gene[i][j]),array.array('H',gene[c][k])))
                matching.append(res)
                message+="<td>"+str(res)+ "</td> \n</tr>\n"
                f.write(str(res)+"\t")
message+= "</tbody>\n</table>\n</details>\n"  

message+= "</html>"  
fh.write(message)
f.close()
fh.close()

