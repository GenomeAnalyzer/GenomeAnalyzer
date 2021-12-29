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
	res = DNA.generating_mRNA(array.array('H', [0,1,0,0,1,1,0,0,0,1,1,0]))
	assert res == "GAUAGC"

# def test_DNA_detecting_genes():
# 	a = DNA.detecting_genes(array.array('I', [0,0,1,1,0,1,0,0,1,1,0,1,1,1,0,0,0,0]))
# 	assert a == [[6,12]]

def test_DNA_generating_amino_acid_chain():
	aa = "NWP"
	# UnicodeDecodeError: 'utf-8' codec can't decode byte 0xf6 in position 4: invalid start byte
	# UnicodeDecodeError: 'utf-8' codec can't decode byte 0xc5 in position 3: invalid continuation byte

	# bin_dna_seq_array = DNA.convert_to_binary(aa,len(aa))

	bin_dna_seq = array.array('H', [0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0])
	# bin_dna_seq = array.array('H', bin_dna_seq_array)
	res = DNA.generating_amino_acid_chain(bin_dna_seq)
	assert aa == res

def test_DNA_detecting_mutations():
	dna_arr = array.array('H', [0,1,0,0,1,1,0,0,0,1,1,0])

	#ar1 = array.array('I', [0 for i in range(3)])
	#ar2 = array.array('I', [0 for i in range(3)])
	#ar3 = array.array('I', [0 for i in range(3)])
	#for i in range(5):
	#	ar1[i] = 0
	#	ar2[i] = 0
	#	ar3[i] = 0
	#res = list([ar1,ar2,ar3])
	#res = list([[0,0,0],[0,0,0],[0,0,0],[0,0,0],[0,0,0]])

	#res = DNA.detecting_mutations(dna_arr)
	res = DNA.detecting_mutations(array.array('H', [0,1,0,1,1,0,0,0,0,0,0,0]))
	assert res[0] == [5,0,5]
	assert res[1] == [0,0,0]
	assert res[2] == [0,0,0]
	assert res[3] == [0,0,0]
	assert res[4] == [0,0,0]

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

