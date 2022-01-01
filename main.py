#coding: utf8
import DNA
import array
import glob
import moduleDNA as m
import os


read = []
sequence = [[]]
mrna = [[[]]]
i = 0
gene = list((()))

chain = [[[]]]

mutation = [[[]]]

matching = [[]]
#gene = [[]]

try:
    # Create target Directory
    os.makedirs("output/sequences")

    print("Output directory created ") 
except FileExistsError:
    print("Output directory already exists")

f = open("output/rapport.txt", "w")
fh = open('output/rapport.html','w')
start = """<html>
<head><style> 
th, td {
  font-size: 10px; 
}
.title {
  font-size: 15px; 
}
/*Style du tableau*/
table, th, td {
  border: 1px solid black;
  border-collapse: collapse;
  border-style: dashed;
}
.title {
  border-style : dashed dashed dashed solid;
  padding-left: 1%;
}
table {
  width:90%;
  margin-left : 5%;
}

details > summary {
  padding: 4px;
  width: 200px;
  background-color: #eeeeee;
  border: none;
  box-shadow: 1px 1px 2px #bbbbbb;
  cursor: help;
}
</style>
</head>

"""

fh.write(start)
for file in glob.glob("fastas/*.fasta"):
    message = "<details><summary>"+str(file)+"</summary>"
    fh.write(message+"<a href=\"sequences/"+str(file).replace("fastas/","")+".html\">"+str(file)+"</a></details>")
    
    f.write("\n######################"+file+"######################")
    fhtmp = open("output/sequences/"+str(file).replace("fastas/","")+'.html','w')

    fhtmp.write(start+ "<h1>"+str(file)+"</h1><a href=\"../rapport.html\" target=\"_blank\"><input type=\"button\" value=\"Retour\"></a>")


    message = "<table>\n<tbody>\n<tr>\n<td class = \"title\">Sequence</td>\n<td class = \"title\">MRNA</td>\n<td class = \"title\">Chain</td>\n<td class = \"title\">Mutation</td>\n</tr>\n"
    

    read = m.read_file(file)

    sequence.append(DNA.convert_to_binary(read,len(read)))
    #print(len(sequence[i]))

    gen =DNA.detecting_genes(array.array('L',sequence[i]))
    gene.append(DNA.detecting_genes(array.array('L',sequence[i])))

    if len(gene[i]) > 0:
      print("gene 1 = "+str(gene[i][0]))
      print("gene 2 = "+str(gen[0]))

    #f.write("\nMRNA\t\t\t\t\t\t\t\tCHAIN\t\t\t\t\t\t\t\tMUTATION\n")
    for j in range (len(gene[i])-1):
        message += "<tr><td>"+str(sequence[i][gene[i][j][0]:(gene[i][j][1])+1])+ "</td>\n"
        res = (DNA.generating_mRNA(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1])])))
        if res:
            mrna.append(res.decode("cp1252", "replace"))
            #f.write(str(res.decode("cp1252", "replace")))
            message+="<td>"+str(res.decode("cp1252", "replace"))+ "</td>\n"
            #print(res.decode("cp1252", "replace"))
        else:
            message+="<td> none </td>\n"
          

        #print('aminoacid ='+str(sequence[i][gene[i][j][0]:(gene[i][j][1])]))
        res2 =  DNA.generating_amino_acid_chain(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1])]))

        if res2:
            chain.append(res2.decode("cp1252", "replace"))
            #f.write("\t\t\t\t\t\t\t\t"+str(res2.decode("cp1252", "replace")))
            message+="<td>"+str(res2.decode("cp1252", "replace"))+ "</td>\n"
        else:
            message+="<td> none </td>\n"            
        #print('ok')
        mutres = DNA.detecting_mutations(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1])]))
        mutation.append(mutres)
        if len(mutres) == 0:
          message+="<td>none</td>\n</tr>\n"
        else:
          message+="<td>"+str(mutres)+ "</td>\n</tr>\n"
        #f.write("\t\t\t\t\t\t\t\t"+str(mutres)+"\n")
    message+= "</tbody>\n</table>\n"

    #print(str(gene))



    f.write("\nMRNA:")
    f.write(str(mrna))
    f.write("\nCHAIN:")
    f.write(str(chain))    
    f.write("\nMUTATION:")
    f.write(str(mutation))
    message += "<table>\n<tr>\n<th class = \"title\">Sequence</th>\n<th class = \"title\">Matching</th>\n</tr>\n<tbody>\n"

    f.write("\nMatching:")
    if (len(gene[i])-1 > 2):
        for j in range((int((len(gene[i]))/2)), -1, -1):
            for k in range((int((len(gene[i]))/2))):     
                message+="<tr><td>"+str(sequence[i][gene[i][j][0]:gene[i][j][1]])+" , "+str(sequence[i][gene[i][k][0]:gene[i][k][1]])+ "</td>\n"
                f.write("["+str(gene[i][j][0])+","+str(gene[i][k][0])+"]=")
                res = (DNA.calculating_matching_score(array.array('H',sequence[i][gene[i][j][0]:gene[i][j][1]]),array.array('H',sequence[i][gene[i][k][0]:gene[i][k][1]])))
                matching.append(res)
                message+="<td>"+str(res)+ "</td> \n</tr>\n"
                f.write(str(res)+"\t")
    message+= "</tbody>\n</table>\n</details>\n </html>"  


    fhtmp.write(message)
    fhtmp.close()
    i=i+1
    if i ==4:
      break


message= " <a href=\"comp.html\">Comparaison séquences</a></details></html>"  
fh.write(message)
fh.close()

fh = open('output/comp.html','w')
fh.write(start+ "<h1>Comparaison entre séquences</h1><a href=\"rapport.html\" target=\"_blank\"><input type=\"button\" value=\"Retour\"></a>")
message = "<table>\n<tr>\n<th class = \"title\">Sequence</th>\n<th class = \"title\">Matching</th>\n</tr>\n<tbody>\n"
f.write("\nMatching:")
for i in range((int((len(gene))/2)), -1, -1):
    for c in range((int((len(gene))/2))):     
        for j in range(int((len(gene[i])))):
            for k in range(int((len(gene[c])))):
                message+="<tr><td>"+str(sequence[i][gene[i][j][0]:gene[i][j][1]])+" , "+str(sequence[c][gene[c][k][0]:gene[c][k][1]])+ "</td>\n"
                f.write("["+str(gene[i][j][0])+","+str(gene[c][k])+"]=")
                res = (DNA.calculating_matching_score(array.array('H',sequence[i][gene[i][j][0]:gene[i][j][1]]),array.array('H',sequence[c][gene[c][k][0]:gene[c][k][1]])))
                matching.append(res)
                message+="<td>"+str(res)+ "</td> \n</tr>\n"
                f.write(str(res)+"\t")
message+= "</tbody>\n</table>\n\n"  

message+= "</html>"  
fh.write(message)
f.close()
fh.close()

