def read_fasta_txt(name):
    gene = ""
    with open(name, "r") as file:
        rows = file.readlines()
        for row in rows:
            gene = gene + row
    return gene

def read_fasta(name):
    gene = ""
    with open(name, "r") as file:
        rows = file.readlines()
        for row in rows[1:]:
            gene = gene + row.strip()
    return gene