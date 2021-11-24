import pytest
import DNA
import array
import moduleDNA as m

###### Verify if C extension in Python is working #####

def test_DNA_convert_to_binary():
	gene = "AATTG"
	print("uh")
	print(type(len(gene)), type(gene))
	res = DNA.compress(len(gene), gene)
	print("uh")
	# printf(res, [0,0,0,0,1,1,1,1,0,1])
	assert res == [0,0,0,0,1,1,1,1,0,1]


def test_DNA_module_read_file():
	res = m.read_file("fastas/test.fasta")
	assert res == "AATTG"

def test_DNA_generating_mRNA():
	res = DNA.generating_mRNA(array.array('i', [0,1,0,0,1,1,0,0,0,1,1,0]))
	assert res == "GAUAGC"

def test_DNA_detecting_genes():
	a = DNA.detecting_genes(array.array('I', [0,0,1,1,1,0,1,1,0,0,0,0]))
	assert a == [[0,6]]

def test_DNA_generating_amino_acid_chain():
	# gene = m.readfile("NC_045512.2.fasta")
	# gene = "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT"
	# gene_bin = DNA.convert_to_binary(gene, len(gene))
	gene_bin = [0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1]
	
	# res = DNA.generating_amino_acid_chain("AUGAGC", 6, 1)
	# assert res == None
	assert None == None

def test_DNA_detecting_mutations():
	res = DNA.detecting_mutations(array.array('I', [0,1,0,0,1,1,0,0,0,1,1,0]))
	assert res

def test_DNA_calculating_matching_score():
	res = DNA.calculating_matching_score(array.array('i', [0,1,0,0,1,1,0,0,0,1,1,0]), array.array('i', [1,0,1,1,0,0,1,1,1,0,0,1]))
	assert res == -1.0

def test_DNA_hamming():
	res = DNA.hamming(0b010011000110, 0b101100111001)
	assert res == 12

if __name__ == "__main__":
	test_DNA_convert_to_binary()
	# test_DNA_generating_amino_acid_chain()