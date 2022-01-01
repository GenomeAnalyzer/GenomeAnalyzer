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
	assert 1 == DNA_bin.get_binary_value(array.array('I', [1]),0)
	# assert leading 0
	assert 0 == DNA_bin.get_binary_value(array.array('I', [1]),17)

  	# Pour chaque binaire de 00000 à 11111, vérifier chaque bit.
	for i in range(32):
		seq_bin = array.array('I', [i])
		assert i%2 == DNA_bin.get_binary_value(seq_bin,0)
		assert int(i/2%2) == DNA_bin.get_binary_value(seq_bin,1)
		assert int(i/4%2) == DNA_bin.get_binary_value(seq_bin,2)
		assert int(i/8%2) == DNA_bin.get_binary_value(seq_bin,3)
		assert int(i/16%2) == DNA_bin.get_binary_value(seq_bin,4)

def test_change_binary_value():
	# Test if the algorithm is OK

  	# 85 = 1010101
	seq_bin = array.array('I', [85])
	# Invert values of seq_bin
	for i in range(0,7): DNA_bin.change_binary_value(seq_bin, i, i%2)
	# 42 = 0101010
	assert 42 == seq_bin[0]
	# Invert back values of seq_bin
	for i in range(0,7): DNA_bin.change_binary_value(seq_bin, i, int((i+1)%2))
	assert 85 == seq_bin[0]

def test_set_binary_array():
	# Test if the algorithm is OK

	assert [85] == DNA_bin.set_binary_array("CCCC", 4)
	assert [170] == DNA_bin.set_binary_array("GGGG", 4)

	assert [255] == DNA_bin.set_binary_array("TTTT", 4)
	assert [0] == DNA_bin.set_binary_array("AAAA", 4)
	assert [0] == DNA_bin.set_binary_array("AAAA", 4)

	seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGA"

	resbin = DNA_bin.set_binary_array(seq_char, len(seq_char))
	assert resbin == [2101911378, 172292753, 4029142153]

def test_xor_binary_array():
  	# Test if the algorithm is OK
	# 00000 xor 11111
	assert [0] == DNA_bin.xor_binary_array(array.array('I', [31]), array.array('I', [31]))
	# 11111 xor 11111
	assert [31] == DNA_bin.xor_binary_array(array.array('I', [0]), array.array('I', [31]))
	# 1010101 xor 0101010
	# res = DNA_bin.xor_binary_array(array.array('I', [42]), array.array('I', [85]))
	# print(res)

def test_popcount_binary_array():
	# Test for all binaries from 00000 to 11111
	for i in range(0,32):
		popc_expected_result = int(i%2) + int(i/2%2) + int(i/4%2) + int(i/8%2) + int(i/16%2)
		popc_result = DNA_bin.popcount_binary_array(array.array('I', [i]))
		assert popc_expected_result == popc_result

	# 9350764 = 001101100111010101110001
	assert 13 == DNA_bin.popcount_binary_array(array.array('I', [9350764]))

def test_mask_binary_array():
	assert 10 == DNA_bin.mask_binary_array(203, 1, 4)


def test_convert_to_binary():
	# Test if the algorithm is OK

	assert [85] == DNA_bin.convert_to_binary("CCCC", 4)
	assert [170] == DNA_bin.convert_to_binary("GGGG", 4)

	assert [255] == DNA_bin.convert_to_binary("TTTT", 4)
	assert [0] == DNA_bin.convert_to_binary("AAAA", 4)
	assert [0] == DNA_bin.convert_to_binary("AAAA", 4)

	seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGA"

	resbin = DNA_bin.convert_to_binary(seq_char, len(seq_char))
	assert resbin == [2101911378, 172292753, 4029142153]

def test_generating_mRNA():
	# Test if the algorithm is OK

	# assert empty
	assert b'' == DNA_bin.generating_mRNA(array.array('I', [0]))
	# 9350764 = 001101100111010101110001 - array size : 1
	assert b'AUGCGUGGGUAG' == DNA_bin.generating_mRNA(array.array('I', [9350764]))
	# 913666358 = 00110110011101010111000100110110, 30065 = 0111010101110001 - array size : 2
	assert b'AUGCGUGGGUAGAUGCGUGGGUAG' == DNA_bin.generating_mRNA(array.array('I', [1821290092, 18263]))
	# assert TypeError for not supported format entries
	with pytest.raises(TypeError):
		DNA_bin.generating_mRNA(None) # no entry
		DNA_bin.generating_mRNA(array.array('H', [12])) # array format double

def test_binary_to_dna():
	assert 0 == 0
	# Test binary to aa conversions

	# --- Test all conversion
	assert b'ATCG' == DNA_bin.binary_to_dna(array.array('I', [156]))
	assert b'ATGCGTGGGTAG' == DNA_bin.binary_to_dna(array.array('I', [9350764]))
	assert b'ATGCGTGGGTAGATGCGTGGGTAG' == DNA_bin.binary_to_dna(array.array('I', [1821290092, 18263]))

	# Test whether the function correctly detects errors:
	with pytest.raises(TypeError):
		DNA_bin.binary_to_dna(None) # no entry
		DNA_bin.binary_to_dna(array.array('H', [12])) # array format double

def test_calculating_matching_score():
	# Test if the algorithm is OK
	# --- With same size
	#  GACCCGAC = 0100101010010010 = {18770}
	#  GGCCAGGC = 0101101000010110 = {26714}
	assert 81.25 == DNA_bin.calculating_matching_score(array.array('I', [18770]), array.array('I', [26714]))
	# --- With different size
	#  GACCCGAC = 0100101010010010 = {18770}
	#  TTTCAGGCTC = 11111110000101101110 = {485503}
	assert 25.0 == DNA_bin.calculating_matching_score(array.array('I', [18770]), array.array('I', [485503]))
	#  TTTCAGGCTT = 11111110000101101111 = {1009791}
	#  GACCTTCGA = 1001010111110010 = {40786}
	assert 35.0 == DNA_bin.calculating_matching_score(array.array('I', [1009791]), array.array('I', [40786]))

	# Test whether the function correctly detects errors:
	# --- NULL error
	with pytest.raises(TypeError):
		DNA_bin.calculating_matching_score(None) # no entry
		DNA_bin.calculating_matching_score(array.array('H', [12])) # array format double

def test_generating_amino_acid_chain():
	assert 0 == 0
	# # Test if the algorithm is OK
	assert b'KKNNR' == DNA_bin.generating_amino_acid_chain(array.array('I', [146868224]))
	assert b'RSSTT' == DNA_bin.generating_amino_acid_chain(array.array('I', [605259304]))
	assert b'TTIMI' == DNA_bin.generating_amino_acid_chain(array.array('I', [481348884]))
	assert b'IEEDD' == DNA_bin.generating_amino_acid_chain(array.array('I', [843718844]))
	assert b'GGGGA' == DNA_bin.generating_amino_acid_chain(array.array('I', [115976842]))
	assert b'AAAVV' == DNA_bin.generating_amino_acid_chain(array.array('I', [775644582]))
	assert b'VVQQH' == DNA_bin.generating_amino_acid_chain(array.array('I', [293871518]))
	assert b'HRRRR' == DNA_bin.generating_amino_acid_chain(array.array('I', [963023473]))
	assert b'PPPPL' == DNA_bin.generating_amino_acid_chain(array.array('I', [232085829]))
	assert b'LLLOO' == DNA_bin.generating_amino_acid_chain(array.array('I', [588240749]))
	assert b'YYOWC' == DNA_bin.generating_amino_acid_chain(array.array('I', [464305363]))
	assert b'CSSSS' == DNA_bin.generating_amino_acid_chain(array.array('I', [928936443]))
	assert b'LLFF' == DNA_bin.generating_amino_acid_chain(array.array('I', [16645071]))
	assert b'MRGOMRGO' == DNA_bin.generating_amino_acid_chain(array.array('I', [1821290092, 18263]))
	# --- Test all the amino acid
	# assert b'KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLLOOYYOWCCSSSSLLFF' == DNA_bin.generating_amino_acid_chain(array.array('I', [146868224, 605259304, 481348884, 843718844, 115976842, 775644582, 293871518, 963023473, 232085829, 588240749, 464305363, 928936443, 16645071]))
	# assert b'KNKNTTTTRSRSIIIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYSSSSCWCLFLFMOOO' == DNA_bin.generating_amino_acid_chain(array.array('I', [79823872, -2096862942,
	# 	-1577991368, 547545866, -1792699787, -1126245655, 1210084514, -752012202, 1001024414, -106443080, -1380064261, -1612777443, 189184]))


	# Test whether the function correctly detects errors:
	with pytest.raises(TypeError):
		DNA_bin.generating_amino_acid_chain(None) # no entry
		DNA_bin.generating_amino_acid_chain(array.array('H', [12])) # array format double

def test_detecting_genes():


   #Test if the algorithm is OK in a basic case: xxxxAUGxxxxUAAxxx
   #The algorithm should detect one gene from the start codon to the stop codon

	assert 6 == DNA_bin.detecting_genes(array.array('i',[963808024, 42]))[0][0]
	assert 28 == DNA_bin.detecting_genes(array.array('i',[963808024, 42]))[0][1]
	assert 1 == len(DNA_bin.detecting_genes(array.array('i',[963808024, 42])))

  	#Test if the algorithm is OK in a non presence of "start/stop" case: xxxxxxxxx
  	#The algorithm should not detect any genes
	assert 0 == len(DNA_bin.detecting_genes(array.array('i',[22369621])))

	
	
	#Test if the algorithm is OK in a multiple "start" case: xxxxAUGxxxxAUGxxxUAAxxx
  	#The algorithm should detect one gene from the 2nd start codon to the stop codon
	assert 30 == DNA_bin.detecting_genes(array.array('i',[732875499, -2036213923]))[0][0]
	assert 48 == DNA_bin.detecting_genes(array.array('i',[732875499, -2036213923]))[0][1]
	assert 1 == len(DNA_bin.detecting_genes(array.array('i',[732875499, -2036213923])))
	
  	#Test if the algorithm is OK in a multiple "stop" case: xxxxAUGxxxxUAAxxxUAAxxx
  	#The algorithm should detect one gene from the start codon to the first stop codon	
	assert 10 == DNA_bin.detecting_genes(array.array('i',[250327787, -2022340747]))[0][0]
	assert 26 == DNA_bin.detecting_genes(array.array('i',[250327787, -2022340747]))[0][1]
	assert 1 == len(DNA_bin.detecting_genes(array.array('i',[250327787, -2022340747])))
	
  	#Test if the algorithm is OK in a multiple gene case: xxxxAUGxxxxUAGxxxAUGxxxUAAxxx
  	#The algorithm should detect two genes
	assert 6 == DNA_bin.detecting_genes(array.array('i',[-469763265, -1612578969, -268435456]))[0][0]
	assert 24 == DNA_bin.detecting_genes(array.array('i',[-469763265, -1612578969, -268435456]))[0][1]
	assert 36 == DNA_bin.detecting_genes(array.array('i',[-469763265, -1612578969, -268435456]))[1][0]
	assert 48 == DNA_bin.detecting_genes(array.array('i',[-469763265, -1612578969, -268435456]))[1][1]
	assert 2 == len(DNA_bin.detecting_genes(array.array('i',[-469763265, -1612578969, -268435456])))

	
def test_detecting_mutations():
  	#GGGTTGCGCGCGTTAAAGGTTTGAAAGGTG = {261725162, 97523700}
  	#Test if sequence 10 to 23 is a mutation zone and no other mutation zone
	assert 13 == DNA_bin.detecting_mutations(array.array('i',[261725162, 97523700]))[0][0]
	assert 10 == DNA_bin.detecting_mutations(array.array('i',[261725162, 97523700]))[0][1]
	assert 23 == DNA_bin.detecting_mutations(array.array('i',[261725162, 97523700]))[0][2]
  	#No other mutations, should not be updated
	assert 0 == DNA_bin.detecting_mutations(array.array('i',[261725162, 97523700]))[1][0]
	assert 0 == DNA_bin.detecting_mutations(array.array('i',[261725162, 97523700]))[1][1]
	assert 0 == DNA_bin.detecting_mutations(array.array('i',[261725162, 97523700]))[1][2]	  

  	#GTTTTGCAAACGTTAAAGGTTTGAAAGGTG = {261102590, 97523700}
  	#Test if no mutation in this sequence
  	 #No possible mutation zones detected, should not be updated	  
	assert 0 == DNA_bin.detecting_mutations(array.array('i',[261102590, 97523700]))[0][0]
	assert 0 == DNA_bin.detecting_mutations(array.array('i',[261102590, 97523700]))[0][1]
	assert 0 == DNA_bin.detecting_mutations(array.array('i',[261102590, 97523700]))[0][2]
  	
	#GGGCCGTTCCGCCCATAGGCCCGGCTAAGA = {-983172758, 17224372}
  	#Test with 3 mutation zones in this sequence
  	#First mutation is updated with right values

	assert 11 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[0][0]
	assert 0 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[0][1]
	assert 11 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[0][2]
  	
	#First mutation is updated with right values

	assert 11 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[1][0]
	assert 16 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[1][1]
	assert 27 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[1][2]
  
   #First mutation is updated with right values

	assert 15 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[2][0]
	assert 34 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[2][1]
	assert 49 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[2][2]


   #No other mutations, should not be updated

	assert 0 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[3][0]
	assert 0 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[3][1]
	assert 0 == DNA_bin.detecting_mutations(array.array('i',[-983172758, 17224372]))[3][2]
  

