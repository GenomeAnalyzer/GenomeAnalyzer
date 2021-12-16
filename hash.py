import json
import itertools

text = """        case CASE :
            aa_seq[temp] = SYMBOL;
            temp++;
            break;
"""

fo = open("codons_list.json")
cdl = json.load(fo)
fo.close()

size = 6

strr = ""

for i in itertools.product([0,1],repeat=size):
    hash = 7000000
    for j in range(size-1,-1,-1):
        mult = 1
        for k in range(size-1,j,-1):
            mult *= 10
        hash += mult * i[j]
    # print(hash)

    shash = str(hash)[1:]
    s2 = "123"
    for j in range(3):
        if shash[2*j:2*j+2] == "00":
            s2 = s2.replace(str(j+1),"A")
        elif shash[2*j:2*j+2] == "11":
            s2 = s2.replace(str(j+1),"T")
        elif shash[2*j:2*j+2] == "10":
            s2 = s2.replace(str(j+1),"C")
        elif shash[2*j:2*j+2] == "01":
            s2 = s2.replace(str(j+1),"G")
    shash = "2" + shash

    # print(shash, s2)

    for item in cdl:
        if item['codon'] == s2.replace("T","U"):
            symbol = item['symbol']
            # print(shash, s2, symbol)

            strr += text.replace("CASE", shash).replace("SYMBOL", f"'{symbol}'")

fw = open("hash_write.txt","w")
fw.write(strr)
fw.close()