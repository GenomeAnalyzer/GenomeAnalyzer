#coding: utf8
import DNA
import array
import glob
import moduleDNA as m
import os
import sys

def main():

  read = list()
  sequence = list()
  i = 0
  gene = list()



  try:
      # Create target Directory
      os.makedirs("output/sequences")

      print("Output directory created ") 
  except FileExistsError:
      print("Output directory already exists")

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
  try:  
    if len(sys.argv) == 1:
      fin = 25698  
    else:    
      if not sys.argv[1].isnumeric():
        raise NameError('nan')
      else:
        fin = sys.argv[1]
  except NameError:
    print("Arg is not a number")
    raise
  fhc = open('output/comp.html','w')

  messagematch = start+ "<h1>Comparaison entre séquences</h1><a href=\"rapport.html\" target=\"_blank\"><input type=\"button\" value=\"Retour\"></a>"

  fh.write(start)
  for file in glob.glob("fastas/*.fasta"):
      
      fhtmp = open("output/sequences/"+str(file).replace("fastas/","")+'.html','w')

      fhtmp.write(start+ "<h1>"+str(file.replace("fastas/","").replace(".fasta",""))+"</h1><a href=\"../rapport.html\" target=\"_blank\"><input type=\"button\" value=\"Retour\"></a>")


      message = "<table>\n<tbody>\n<tr>\n<td class = \"title\">Sequence</td>\n<td class = \"title\">MRNA</td>\n<td class = \"title\">Chain</td>\n<td class = \"title\">Mutation</td>\n</tr>\n"
      

      read = m.read_file(file)


      sequence.append(DNA.convert_to_binary(read,len(read)))

      gene.append(DNA.detecting_genes(array.array('L',sequence[i])))


      fh.write("<details><summary>"+str(file.replace("fastas/","").replace(".fasta",""))+"</summary><a href=\"sequences/"+str(file).replace("fastas/","")+".html\">"+str(file.replace("fastas/","").replace(".fasta",""))+"- Genes: "+str(len(gene[i]))+"</a></details>")

      for j in range (len(gene[i])-1):
          message += "<tr><td>"+str(sequence[i][gene[i][j][0]:(gene[i][j][1])+1])+ "</td>\n"
          res = (DNA.generating_mRNA(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1]+1)])))
          if res:
              message+="<td>"+str(res.decode("cp1252", "replace"))+ "</td>\n"
          else:
              message+="<td> none </td>\n"
            

          res2 =  DNA.generating_amino_acid_chain(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1]+1)]))

          if res2:
              message+="<td>"+str(res2.decode("cp1252", "replace"))+ "</td>\n"
          else:
              message+="<td> none </td>\n"            
          mutres = DNA.detecting_mutations(array.array('H',sequence[i][gene[i][j][0]:(gene[i][j][1]+1)]))
          if len(mutres) == 0:
            message+="<td>none</td>\n</tr>\n"
          else:
            message+="<td>"+str(mutres)+ "</td>\n</tr>\n"
      message+= "</tbody>\n</table>\n"


      message += "<table>\n<tr>\n<th class = \"title\">Sequence</th>\n<th class = \"title\">Matching</th>\n</tr>\n<tbody>\n"
      if (len(gene[i])-1 > 2):
          for j in range((int((len(gene[i]))/2)), -1, -1):
              for k in range((int((len(gene[i]))/2))):     
                  message+="<tr><td>"+str(sequence[i][gene[i][j][0]:gene[i][j][1]])+" , "+str(sequence[i][gene[i][k][0]:gene[i][k][1]])+ "</td>\n"
                  res = (DNA.calculating_matching_score(array.array('H',sequence[i][gene[i][j][0]:gene[i][j][1]+1]),array.array('H',sequence[i][gene[i][k][0]:gene[i][k][1]+1])))
                  message+="<td>"+str(res)+ "</td> \n</tr>\n"
      message+= "</tbody>\n</table>\n</details>\n </html>"  


      fhtmp.write(message)
      fhtmp.close()

      for c in range(i-1, -1, -1): 
            fhtmp2 = open("output/sequences/cmp"+str(i)+"-"+str(c)+'.html','w')

            messagematch+="<details><summary>Sequence "+str(i)+" - "+str(c)+"</summary><a href=\"sequences/cmp"+str(i)+"-"+str(c)+".html\">Comparaison "+str(i)+"-"+str(c)+"</a></details>\n"
            msgtmp = "<table>\n<tr>\n<th class = \"title\">Sequence</th>\n<th class = \"title\">Matching</th>\n</tr>\n<tbody>"

            for j in range(int((len(gene[i])))):
                for k in range(int((len(gene[c])))):
                      msgtmp+="<tr><td>"+str(sequence[i][gene[i][j][0]:gene[i][j][1]])+" , "+str(sequence[c][gene[c][k][0]:gene[c][k][1]])+ "</td>\n"
                      res = (DNA.calculating_matching_score(array.array('H',sequence[i][gene[i][j][0]:gene[i][j][1]+1]),array.array('H',sequence[c][gene[c][k][0]:gene[c][k][1]+1])))
                      msgtmp+="<td>"+str(res)+ "</td> \n</tr>\n"
            fhtmp2.write(msgtmp)
            fhtmp2.close()
      i=i+1
      if i == int(fin):
        break


  message= " <a href=\"comp.html\">Comparaison séquences</a></details></html>"  
  fh.write(message)
  fh.close()


  messagematch+= "</tbody>\n</table>\n\n</html>"  

  fhc.write(messagematch)
  fhc.close()


if __name__ == "__main__":
    main()
