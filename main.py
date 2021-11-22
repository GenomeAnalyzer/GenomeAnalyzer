import DNA
import array
import numpy as np
import faulthandler
import sys
import ctypes 
import glob
import moduleDNA as m



for file in glob.glob("fastas/*.fasta"):
    test = []
    #print(str(file))
    test = m.read_file(file)
#print(test)
#    print(str(len(test)))
    bin = DNA.convert_to_binary(test,len(test)*2)
    #ok = np.ctypeslib.as_array(truc, shape=(3072))
#print(bin)

