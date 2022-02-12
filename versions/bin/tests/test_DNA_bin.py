import pytest
import DNA_bin
import array
import moduleDNA as m
import ctypes

int_SIZE = 31

###### Verify if C extension in Python is working #####


def test_get_binary_value():
	# Test if the algorithm is OK
	# 1 = 0000000000000000000000000000001
	assert 1 == DNA_bin.get_binary_value(array.array('l', [1]),0)
	# assert leading 0
	assert 0 == DNA_bin.get_binary_value(array.array('l', [1]),17)

  	# Pour chaque binaire de 00000 à 11111, vérifier chaque bit.
	for i in range(32):
		seq_bin = array.array('l', [i])
		assert i%2 == DNA_bin.get_binary_value(seq_bin,0)
		assert int(i/2%2) == DNA_bin.get_binary_value(seq_bin,1)
		assert int(i/4%2) == DNA_bin.get_binary_value(seq_bin,2)
		assert int(i/8%2) == DNA_bin.get_binary_value(seq_bin,3)
		assert int(i/16%2) == DNA_bin.get_binary_value(seq_bin,4)

def test_convert_to_binary():
	# Test if the algorithm is OK

	assert [170] == DNA_bin.convert_to_binary("CCCC", 4, 1)
	assert [85] == DNA_bin.convert_to_binary("GGGG", 4, 1)

	assert [255] == DNA_bin.convert_to_binary("TTTT", 4, 1)
	assert [0] == DNA_bin.convert_to_binary("AAAA", 4, 1)

	seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCT\nTCGAGACCTTCGA"

	resbin = DNA_bin.convert_to_binary(seq_char, len(seq_char), 2)
	assert resbin == [5402369836413063467, 59845604]

def test_generating_mRNA():
	# Test if the algorithm is OK

	# assert empty
	assert b'' == DNA_bin.generating_mRNA(array.array('l', [0]), 0, 0)
	# 9350764 = 001101100111010101110001 - array size : 1
	assert b'AUGCGUGGGUAG' == DNA_bin.generating_mRNA(array.array('l', [9350764]), 0, 24)
	# assert TypeError for not supported format entries
	with pytest.raises(TypeError):
		DNA_bin.generating_mRNA(None) # no entry
		DNA_bin.generating_mRNA(array.array('H', [12])) # array format double

def test_binary_to_dna():
	assert 0 == 0
	# Test binary to aa conversions

	# --- Test all conversion
	assert b'ATCG' == DNA_bin.binary_to_dna(array.array('l', [156]))
	assert b'ATGCGTGGGTAG' == DNA_bin.binary_to_dna(array.array('l', [9350764]))
	assert b'ATGCGTGGGTAGATGCAAAAAAAA' == DNA_bin.binary_to_dna(array.array('l', [1821290092, 18263]))

	# Test whether the function correctly detects errors:
	with pytest.raises(TypeError):
		DNA_bin.binary_to_dna(None) # no entry
		DNA_bin.binary_to_dna(array.array('H', [12])) # array format double

def test_calculating_matching_score():
	# Test if the algorithm is OK
	# --- With same size
	#  GACCCGAC = 0100101010010010 = {18770}
	#  GGCCAGGC = 0101101000010110 = {26714}
	assert 81.25 == DNA_bin.calculating_matching_score(array.array('l', [18770]), 0, 16, array.array('l', [26714]), 0, 16)
	# --- With different size
	#  GACCCGAC = 0100101010010010 = {18770}
	#  TTTCAGGCTC = 11111110000101101110 = {485503}
	assert 55.0 == DNA_bin.calculating_matching_score(array.array('l', [18770]), 0, 16, array.array('l', [485503]), 0, 20)
	#  TTTCAGGCTT = 11111110000101101111 = {1009791}
	#  GACCTTCGA = 1001010111110010 = {40786}
	assert 25.0 == DNA_bin.calculating_matching_score(array.array('l', [1009791]), 0, 20, array.array('l', [40786]), 0, 16)

	# Test whether the function correctly detects errors:
	# --- NULL error
	with pytest.raises(TypeError):
		DNA_bin.calculating_matching_score(None) # no entry
		DNA_bin.calculating_matching_score(array.array('H', [12])) # array format double

def test_generating_amino_acid_chain():
	assert 0 == 0
	# # Test if the algorithm is OK
	assert b'KKNNRKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [146868224]), 0, 60)
	assert b'RSSTTKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [605259304]), 0, 60)
	assert b'TTIMIKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [481348884]), 0, 60)
	assert b'IEEDDKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [843718844]), 0, 60)
	assert b'GGGGAKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [115976842]), 0, 60)
	assert b'AAAVVKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [775644582]), 0, 60)
	assert b'VVQQHKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [293871518]), 0, 60)
	assert b'HRRRRKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [963023473]), 0, 60)
	assert b'PPPPLKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [232085829]), 0, 60)
	assert b'LLLOOKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [588240749]), 0, 60)
	assert b'YYOWCKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [464305363]), 0, 60)
	assert b'CSSSSKKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [928936443]), 0, 60)
	assert b'LLFF' == DNA_bin.generating_amino_acid_chain(array.array('l', [16645071]), 0, 24)
	assert b'MRGOMQKKKK' == DNA_bin.generating_amino_acid_chain(array.array('l', [1821290092, 18263]), 0, 60)


	# Test whether the function correctly detects errors:
	with pytest.raises(TypeError):
		DNA_bin.generating_amino_acid_chain(None) # no entry
		DNA_bin.generating_amino_acid_chain(array.array('H', [12])) # array format double

def test_detecting_genes():

   #Test if the algorithm is OK in a basic case: xxxxAUGxxxxUAAxxx
   #The algorithm should detect one gene from the start codon to the stop codon

	assert 6 == DNA_bin.detecting_genes(array.array('l',[963808024, 42]))[0][0]
	assert 33 == DNA_bin.detecting_genes(array.array('l',[963808024, 42]))[0][1]
	assert 1 == len(DNA_bin.detecting_genes(array.array('l',[963808024, 42])))

  	#Test if the algorithm is OK in a non presence of "start/stop" case: xxxxxxxxx
  	#The algorithm should not detect any genes
	assert 0 == len(DNA_bin.detecting_genes(array.array('l',[22369621])))

	
	
	#Test if the algorithm is OK in a multiple "start" case: xxxxAUGxxxxAUGxxxUAAxxx
  	#The algorithm should detect one gene from the 2nd start codon to the stop codon
	assert 62 == DNA_bin.detecting_genes(array.array('l',[732875499, -2036213923]))[0][0]
	assert 85 == DNA_bin.detecting_genes(array.array('l',[732875499, -2036213923]))[0][1]
	assert 1 == len(DNA_bin.detecting_genes(array.array('l',[732875499, -2036213923])))
	
  	#Test if the algorithm is OK in a multiple "stop" case: xxxxAUGxxxxUAAxxxUAAxxx
  	#The algorithm should detect one gene from the start codon to the first stop codon	
	assert 10 == DNA_bin.detecting_genes(array.array('l',[250327787, -2022340747]))[0][0]
	assert 31 == DNA_bin.detecting_genes(array.array('l',[250327787, -2022340747]))[0][1]
	assert 2 == len(DNA_bin.detecting_genes(array.array('l',[250327787, -2022340747])))
	
  	#Test if the algorithm is OK in a multiple gene case: xxxxAUGxxxxUAGxxxAUGxxxUAAxxx
  	#The algorithm should detect two genes
	assert 6 == DNA_bin.detecting_genes(array.array('l',[-469763265, -1612578969, -268435456]))[0][0]
	assert 29 == DNA_bin.detecting_genes(array.array('l',[-469763265, -1612578969, -268435456]))[0][1]
	assert 68 == DNA_bin.detecting_genes(array.array('l',[-469763265, -1612578969, -268435456]))[1][0]
	assert 85 == DNA_bin.detecting_genes(array.array('l',[-469763265, -1612578969, -268435456]))[1][1]
	assert 2 == len(DNA_bin.detecting_genes(array.array('l',[-469763265, -1612578969, -268435456])))

	
def test_detecting_mutations():

  	# GGGTTGCGCGCGTTAAAGGTTTGAAAGGTG = {261725162, 97523700}
  	# Test if sequence 10 to 23 is a mutation zone and no other mutation zone
	assert 13 == DNA_bin.detecting_mutations(array.array('l',[261725162, 97523700]),0,60)[0][0]
	assert 10 == DNA_bin.detecting_mutations(array.array('l',[261725162, 97523700]),0,60)[0][1]
	assert 24 == DNA_bin.detecting_mutations(array.array('l',[261725162, 97523700]),0,60)[0][2]
	
	#GGGCCGTTCCGCCCATAGGCCCGGCTAAGA = {-983172758, 17224372}
  	#Test with 3 mutation zones in this sequence
  	#First mutation is updated with right values
	assert 11 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]),0,60)[0][0]
	assert 0 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]),0,60)[0][1]
	assert 12 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]),0,60)[0][2]

	#Second mutation is updated with right values
	assert 11 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]),0,60)[1][0]
	assert 16 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]),0,60)[1][1]
	assert 28 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]),0,60)[1][2]

