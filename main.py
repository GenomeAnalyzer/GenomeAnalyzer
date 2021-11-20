import DNA
import array
import numpy as np
import faulthandler
import sys
import ctypes 

def read_file(name):
    f = open(name,"r")
    f.readline()
    gene = ""
    while 1:
        char = f.read(1)
        if not char:
            f.close()
            break
        #print(char)
        gene= gene +char
    f.close()
    return gene

test = []
test = read_file("fastas/6X2G_E.fasta")

print(test)
faulthandler.enable()
truc = DNA.convert_to_binary(test+'\n')
#ok = np.ctypeslib.as_array(truc, shape=(3072))
print(truc)
#print("Ok final : "+ ok+"\n")

