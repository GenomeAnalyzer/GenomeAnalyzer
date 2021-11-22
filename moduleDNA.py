def read_file(name):
    f = open(name,"r")
    f.readline()
    gene = ""
    while 1:
        char = f.readline()
        if not char:
            f.close()
            break
        
        #print(ascii(char.replace(" ","").replace("\n","")))
        gene= gene +char.replace(" ","").replace("\n","")
    f.close()
    return gene