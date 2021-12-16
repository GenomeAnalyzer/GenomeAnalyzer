import pytest
import DNA
import array
import moduleDNA as m

###### Verify if C extension in Python is working #####

def test_DNA_convert_to_binary():
	res = DNA.convert_to_binary("AATTG",5)
	assert res == [0,0,0,0,1,1,1,1,0,1]

def test_DNA_module_read_file():
	res = m.read_file("fastas/test.fasta")
	assert res == "AATTG"

def test_DNA_generating_mRNA():
	res = DNA.generating_mRNA(array.array('i', [0,1,0,0,1,1,0,0,0,1,1,0]))
	assert res == "GAUAGC"

def test_DNA_detecting_genes():
	a = DNA.detecting_genes(array.array('I', [0,0,1,1,0,1,0,0,1,1,0,1,1,1,0,0,0,0]))
	assert a == [[6,12]]

#Maybe not good value, I don't know
def test_DNA_generating_amino_acid_chain():
	res = DNA.generating_amino_acid_chain("AUGAGC", 6, 1)
	assert res == None

def test_DNA_detecting_mutations():
	res = DNA.detecting_mutations(array.array('I', [0,1,0,0,1,1,0,0,0,1,1,0]))
	assert res

def test_DNA_calculating_matching_score():
	size = 20
	ar1 = array.array('H', [0 for i in range(size)])
	ar2 = array.array('H', [1 for i in range(size)])
	for i in range(size):
		ar1[i] = 1
		assert (i+1)*100/size == DNA.calculating_matching_score(ar1, ar2)
	# Warning here, (i+1)*100/size might generate a floating number with low precision error, which might result in a false assert.

def test_DNA_hamming():
	res = DNA.hamming(0b010011000110, 0b101100111001)
	assert res == 12

