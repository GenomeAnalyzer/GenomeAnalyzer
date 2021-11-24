def read_file(name):
    f = open(name,"r")
    f.readline()
    gene = ""
    while 1:
        char = f.readline()
        if not char:
            f.close()
            break
        gene= gene +char.replace(" ","").replace("\n","")
    f.close()
    return gene

def readfile(name):
    fo = open(name, "r")
    foo = fo.readlines()
    fo.close()
    gene = "".join([line[:-1].replace(" ","") for line in foo[1:]])
    return gene

if __name__ == "__main__":
    gene = readfile("NC_045512.2.fasta")
    print(gene)