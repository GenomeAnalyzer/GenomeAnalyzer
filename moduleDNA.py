def read_file(name):
    f = open(name,"r")
    f.readline()
    gene = ""
    while 1:
        char = f.readline()
        if not char:
            f.close()
            break
        gene= gene +char.strip()
    f.close()
    return gene

def get_sequence(sequence,start,end):
    print(str(end-start)) 
    print(str(len(sequence)))

    real_start = list()
    real_end = list()
    mid = list()

    gene = list()

    nb = len(sequence[int((start/31)):int(end/31)+1])


    start_case = sequence[int(start/31)]
    l = list()
    print(str(bin(start_case.zfill(31))))
    for x in bin(start_case.zfill(31))[2:]:
        if x.isdigit():
            l.append(int(x))
        
    l.reverse()
    print("reverse ="+str(l) )
    #Si le début n'est pas un multiple de 31 alors c'est dans le millieu d'une case de la liste
    if start%31 != 0 :

        gene.append(l[((start)%31):])

    else:
        gene.append(l)



    if nb > 1:
        end_case = sequence[int(end/31)]
        print(str(bin(end_case.zfill(31))))

        l = list()
        for x in bin(end_case.zfill(31))[2:]:
            if x.isdigit():
                l.append(int(x))
                
        l.reverse()        
        if end%31 != 0 :
                #print(str(int(end/31)))
            if nb > 2:
                gene.append(sequence[int((start/31)):int(end/31)])
                tmp_l = list()
                for i in range(int(start/31),int(end/31)):
                    for x in bin(sequence[i])[2:]:
                        if x.isdigit():
                            tmp_l.append(int(x))
                tmp_l.reverse()
                print("tmp = "+str(tmp_l))
                gene.append(l[:((end+30)%31)])
            else:
                gene.append(l[:((end+30)%31)])
        else:
            if nb > 2:
                gene.append(sequence[int((start/31)):int(end/31)])
                gene.append(l)

            else:
                gene.append(l)
           
    else:
        
        mid = [0]
        real_end = [0]

    print(str(real_start))
    #reconversion après 
    res_start = 0
    for num in real_start:
        res_start = (res_start << 1) | num

    #print("res start"+str(res_start))
    return gene
