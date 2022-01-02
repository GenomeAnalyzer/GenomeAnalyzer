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

def test_change_binary_value():
	# Test if the algorithm is OK

  	# 85 = 1010101
	seq_bin = array.array('l', [85])
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
	# assert resbin == [2101911378, 172292753, 4029142153] # unsigned int
	assert resbin == [2101911378, 172292753, -265825143] # long int

def test_xor_binary_array():
  	# Test if the algorithm is OK
	# 00000 xor 11111
	assert [0] == DNA_bin.xor_binary_array(array.array('l', [31]), array.array('l', [31]))
	# 11111 xor 11111
	assert [31] == DNA_bin.xor_binary_array(array.array('l', [0]), array.array('l', [31]))
	# 1010101 xor 0101010
	# res = DNA_bin.xor_binary_array(array.array('l', [42]), array.array('l', [85]))
	# print(res)

def test_popcount_binary_array():
	# Test for all binaries from 00000 to 11111
	for i in range(0,32):
		popc_expected_result = int(i%2) + int(i/2%2) + int(i/4%2) + int(i/8%2) + int(i/16%2)
		popc_result = DNA_bin.popcount_binary_array(array.array('l', [i]))
		assert popc_expected_result == popc_result

	# 9350764 = 001101100111010101110001
	assert 13 == DNA_bin.popcount_binary_array(array.array('l', [9350764]))

def test_mask_binary_array():
	assert 5 == DNA_bin.mask_binary_array(203, 1, 4)

	intsize = int_SIZE + 1
	size = 511
	
	arr = []
	i = size
	while i >=0:
		arr.append(i)
		i = i-1

	arr = array.array('l', arr)
	# print(len(arr), arr)
	# input()

	for it in range(0,size):
		pos = size - it
		# print(it, pos, arr[pos])

		# retourner le Xième bit
		assert it%2 == DNA_bin.mask_binary_array(arr[pos], 0,1)
		assert int(it/2)%2 == DNA_bin.mask_binary_array(arr[pos], 1,1)
		assert int(it/4)%2 == DNA_bin.mask_binary_array(arr[pos], 2,1)
		assert int(it/8)%2 == DNA_bin.mask_binary_array(arr[pos], 3,1)
		assert int(it/16)%2 == DNA_bin.mask_binary_array(arr[pos], 4,1)
		# Retourner tout le nombre
		assert it == DNA_bin.mask_binary_array(arr[pos], 0,intsize)
		# retourner plusieurs bitss
		assert int(it/2)%2 * 2 + it%2 == DNA_bin.mask_binary_array(arr[pos], 0,2)
		assert int(it/4)%2 * 4 + int(it/2)%2 * 2 + it%2 == DNA_bin.mask_binary_array(arr[pos], 0,3)
		assert int(it/4)%2 * 2 + int(it/2)%2 == DNA_bin.mask_binary_array(arr[pos], 1,2)
		assert int(it/16)%2 * 4 + int(it/8)%2 * 2 + int(it/4)%2 == DNA_bin.mask_binary_array(arr[pos], 2,3)

def test_get_piece_binary_array():
	intsize = int_SIZE + 1
	size = 31

	arr = array.array('l', [0b11111111111111111111111111111111, 0b101000111000101, 0b101000111000101, 0b101000111000101, 0b101000111000101, 31, 31])
	# arr = array.array('l', [0b101000111000101, 0b101000111000101, 0b101000111000101, 0b101000111000101, 31, 31])
	assert [2, 2147494114, 2147494114] == DNA_bin.get_piece_binary_array(arr, 2*intsize+1, 2*intsize+2)
	assert [1, 20933, 20933] == DNA_bin.get_piece_binary_array(arr, 2*intsize, 2*intsize+2)
	assert [17, 1073747057, 1073741831, 3221225479] == DNA_bin.get_piece_binary_array(arr, 2, 3*intsize+5)

	# Tests récupérations premiers entiers
	assert [15] == DNA_bin.get_piece_binary_array(arr, 0, 4)
	assert [31] == DNA_bin.get_piece_binary_array(arr, 0, intsize)

	assert [0b101000111000101] == DNA_bin.get_piece_binary_array(arr, 2*intsize, intsize)
	assert [0b101000111000101, 31] == DNA_bin.get_piece_binary_array(arr, intsize, 2*intsize)
	assert [0b11111111111111111111111111111111] == DNA_bin.get_piece_binary_array(arr, 6*intsize, intsize)

	# SURCHARGE DES TESTS AVEC LES TESTS DE mask_binary_array
	intsize = int_SIZE + 1
	size = 511
	
	arr = []
	i = size
	while i >=0:
		arr.append(i)
		i = i-1

	arr = array.array('l', arr)
	# print(len(arr), arr)

	for it in range(0,size):
		pos = size - it
		# print("")
		# print("it : ", it, "pos : ", pos, "arr[pos] : ", arr[pos])

		# retourner le Xième bit
		assert [it%2] == DNA_bin.get_piece_binary_array(arr, it*intsize + 0,1)
		assert [int(it/2)%2] == DNA_bin.get_piece_binary_array(arr, it*intsize + 1,1)
		assert [int(it/4)%2] == DNA_bin.get_piece_binary_array(arr, it*intsize + 2,1)
		assert [int(it/8)%2] == DNA_bin.get_piece_binary_array(arr, it*intsize + 3,1)
		assert [int(it/16)%2] == DNA_bin.get_piece_binary_array(arr, it*intsize + 4,1)
		# Retourner tout le nombre
		assert [it] == DNA_bin.get_piece_binary_array(arr, it*intsize + 0,intsize)
		# retourner plusieurs bitss
		assert [int(it/2)%2 * 2 + it%2] == DNA_bin.get_piece_binary_array(arr, it*intsize + 0,2)
		assert [int(it/4)%2 * 4 + int(it/2)%2 * 2 + it%2] == DNA_bin.get_piece_binary_array(arr, it*intsize + 0,3)
		assert [int(it/4)%2 * 2 + int(it/2)%2] == DNA_bin.get_piece_binary_array(arr, it*intsize + 1,2)
		assert [int(it/16)%2 * 4 + int(it/8)%2 * 2 + int(it/4)%2] == DNA_bin.get_piece_binary_array(arr, it*intsize + 2,3)

		if it < size -1 : # de 0 à size -2
			assert [it+1, it] == DNA_bin.get_piece_binary_array(arr, it*intsize,2*intsize)
			assert [(it+2)%2, it+1, it] == DNA_bin.get_piece_binary_array(arr, it*intsize,2*intsize+1)


def test_convert_to_binary():
	# Test if the algorithm is OK

	assert [85] == DNA_bin.convert_to_binary("CCCC", 4)
	assert [170] == DNA_bin.convert_to_binary("GGGG", 4)

	assert [255] == DNA_bin.convert_to_binary("TTTT", 4)
	assert [0] == DNA_bin.convert_to_binary("AAAA", 4)
	assert [0] == DNA_bin.convert_to_binary("AAAA", 4)

	seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGA"

	resbin = DNA_bin.convert_to_binary(seq_char, len(seq_char))
	# assert resbin == [2101911378, 172292753, 4029142153] # if unsigned int
	assert resbin == [2101911378, 172292753, -265825143] # if long int

def test_generating_mRNA():
	# Test if the algorithm is OK

	# assert empty
	assert b'' == DNA_bin.generating_mRNA(array.array('l', [0]))
	# 9350764 = 001101100111010101110001 - array size : 1
	assert b'AUGCGUGGGUAG' == DNA_bin.generating_mRNA(array.array('l', [9350764]))
	# 913666358 = 00110110011101010111000100110110, 30065 = 0111010101110001 - array size : 2
	assert b'AUGCGUGGGUAGAUGCGUGGGUAG' == DNA_bin.generating_mRNA(array.array('l', [1821290092, 18263]))
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
	assert b'ATGCGTGGGTAGATGCGTGGGTAG' == DNA_bin.binary_to_dna(array.array('l', [1821290092, 18263]))

	# Test whether the function correctly detects errors:
	with pytest.raises(TypeError):
		DNA_bin.binary_to_dna(None) # no entry
		DNA_bin.binary_to_dna(array.array('H', [12])) # array format double

def test_calculating_matching_score():
	# Test if the algorithm is OK
	# --- With same size
	#  GACCCGAC = 0100101010010010 = {18770}
	#  GGCCAGGC = 0101101000010110 = {26714}
	assert 81.25 == DNA_bin.calculating_matching_score(array.array('l', [18770]), array.array('l', [26714]))
	# --- With different size
	#  GACCCGAC = 0100101010010010 = {18770}
	#  TTTCAGGCTC = 11111110000101101110 = {485503}
	assert 25.0 == DNA_bin.calculating_matching_score(array.array('l', [18770]), array.array('l', [485503]))
	#  TTTCAGGCTT = 11111110000101101111 = {1009791}
	#  GACCTTCGA = 1001010111110010 = {40786}
	assert 35.0 == DNA_bin.calculating_matching_score(array.array('l', [1009791]), array.array('l', [40786]))

	# Test whether the function correctly detects errors:
	# --- NULL error
	with pytest.raises(TypeError):
		DNA_bin.calculating_matching_score(None) # no entry
		DNA_bin.calculating_matching_score(array.array('H', [12])) # array format double

def test_generating_amino_acid_chain():
	assert 0 == 0
	# # Test if the algorithm is OK
	assert b'KKNNR' == DNA_bin.generating_amino_acid_chain(array.array('l', [146868224]))
	assert b'RSSTT' == DNA_bin.generating_amino_acid_chain(array.array('l', [605259304]))
	assert b'TTIMI' == DNA_bin.generating_amino_acid_chain(array.array('l', [481348884]))
	assert b'IEEDD' == DNA_bin.generating_amino_acid_chain(array.array('l', [843718844]))
	assert b'GGGGA' == DNA_bin.generating_amino_acid_chain(array.array('l', [115976842]))
	assert b'AAAVV' == DNA_bin.generating_amino_acid_chain(array.array('l', [775644582]))
	assert b'VVQQH' == DNA_bin.generating_amino_acid_chain(array.array('l', [293871518]))
	assert b'HRRRR' == DNA_bin.generating_amino_acid_chain(array.array('l', [963023473]))
	assert b'PPPPL' == DNA_bin.generating_amino_acid_chain(array.array('l', [232085829]))
	assert b'LLLOO' == DNA_bin.generating_amino_acid_chain(array.array('l', [588240749]))
	assert b'YYOWC' == DNA_bin.generating_amino_acid_chain(array.array('l', [464305363]))
	assert b'CSSSS' == DNA_bin.generating_amino_acid_chain(array.array('l', [928936443]))
	assert b'LLFF' == DNA_bin.generating_amino_acid_chain(array.array('l', [16645071]))
	assert b'MRGOMRGO' == DNA_bin.generating_amino_acid_chain(array.array('l', [1821290092, 18263]))
	# --- Test all the amino acid
	# assert b'KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLLOOYYOWCCSSSSLLFF' == DNA_bin.generating_amino_acid_chain(array.array('l', [146868224, 605259304, 481348884, 843718844, 115976842, 775644582, 293871518, 963023473, 232085829, 588240749, 464305363, 928936443, 16645071]))
	# assert b'KNKNTTTTRSRSIIIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYSSSSCWCLFLFMOOO' == DNA_bin.generating_amino_acid_chain(array.array('l', [79823872, -2096862942,
	# 	-1577991368, 547545866, -1792699787, -1126245655, 1210084514, -752012202, 1001024414, -106443080, -1380064261, -1612777443, 189184]))


	# Test whether the function correctly detects errors:
	with pytest.raises(TypeError):
		DNA_bin.generating_amino_acid_chain(None) # no entry
		DNA_bin.generating_amino_acid_chain(array.array('H', [12])) # array format double

def test_detecting_genes():


   #Test if the algorithm is OK in a basic case: xxxxAUGxxxxUAAxxx
   #The algorithm should detect one gene from the start codon to the stop codon

	assert 6 == DNA_bin.detecting_genes(array.array('l',[963808024, 42]))[0][0]
	assert 28 == DNA_bin.detecting_genes(array.array('l',[963808024, 42]))[0][1]
	assert 1 == len(DNA_bin.detecting_genes(array.array('l',[963808024, 42])))

  	#Test if the algorithm is OK in a non presence of "start/stop" case: xxxxxxxxx
  	#The algorithm should not detect any genes
	assert 0 == len(DNA_bin.detecting_genes(array.array('l',[22369621])))

	
	
	#Test if the algorithm is OK in a multiple "start" case: xxxxAUGxxxxAUGxxxUAAxxx
  	#The algorithm should detect one gene from the 2nd start codon to the stop codon
	assert 30 == DNA_bin.detecting_genes(array.array('l',[732875499, -2036213923]))[0][0]
	assert 48 == DNA_bin.detecting_genes(array.array('l',[732875499, -2036213923]))[0][1]
	assert 1 == len(DNA_bin.detecting_genes(array.array('l',[732875499, -2036213923])))
	
  	#Test if the algorithm is OK in a multiple "stop" case: xxxxAUGxxxxUAAxxxUAAxxx
  	#The algorithm should detect one gene from the start codon to the first stop codon	
	assert 10 == DNA_bin.detecting_genes(array.array('l',[250327787, -2022340747]))[0][0]
	assert 26 == DNA_bin.detecting_genes(array.array('l',[250327787, -2022340747]))[0][1]
	assert 1 == len(DNA_bin.detecting_genes(array.array('l',[250327787, -2022340747])))
	
  	#Test if the algorithm is OK in a multiple gene case: xxxxAUGxxxxUAGxxxAUGxxxUAAxxx
  	#The algorithm should detect two genes
	assert 6 == DNA_bin.detecting_genes(array.array('l',[-469763265, -1612578969, -268435456]))[0][0]
	assert 24 == DNA_bin.detecting_genes(array.array('l',[-469763265, -1612578969, -268435456]))[0][1]
	assert 36 == DNA_bin.detecting_genes(array.array('l',[-469763265, -1612578969, -268435456]))[1][0]
	assert 48 == DNA_bin.detecting_genes(array.array('l',[-469763265, -1612578969, -268435456]))[1][1]
	assert 2 == len(DNA_bin.detecting_genes(array.array('l',[-469763265, -1612578969, -268435456])))

	
def test_detecting_mutations():
  	#GGGTTGCGCGCGTTAAAGGTTTGAAAGGTG = {261725162, 97523700}
  	#Test if sequence 10 to 23 is a mutation zone and no other mutation zone
	assert 13 == DNA_bin.detecting_mutations(array.array('l',[261725162, 97523700]))[0][0]
	assert 10 == DNA_bin.detecting_mutations(array.array('l',[261725162, 97523700]))[0][1]
	assert 23 == DNA_bin.detecting_mutations(array.array('l',[261725162, 97523700]))[0][2]
  	#No other mutations, should not be updated
	assert 0 == DNA_bin.detecting_mutations(array.array('l',[261725162, 97523700]))[1][0]
	assert 0 == DNA_bin.detecting_mutations(array.array('l',[261725162, 97523700]))[1][1]
	assert 0 == DNA_bin.detecting_mutations(array.array('l',[261725162, 97523700]))[1][2]	  

  	#GTTTTGCAAACGTTAAAGGTTTGAAAGGTG = {261102590, 97523700}
  	#Test if no mutation in this sequence
  	 #No possible mutation zones detected, should not be updated	  
	assert 0 == DNA_bin.detecting_mutations(array.array('l',[261102590, 97523700]))[0][0]
	assert 0 == DNA_bin.detecting_mutations(array.array('l',[261102590, 97523700]))[0][1]
	assert 0 == DNA_bin.detecting_mutations(array.array('l',[261102590, 97523700]))[0][2]
  	
	#GGGCCGTTCCGCCCATAGGCCCGGCTAAGA = {-983172758, 17224372}
  	#Test with 3 mutation zones in this sequence
  	#First mutation is updated with right values

	assert 11 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[0][0]
	assert 0 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[0][1]
	assert 11 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[0][2]
  	
	#First mutation is updated with right values

	assert 11 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[1][0]
	assert 16 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[1][1]
	assert 27 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[1][2]
  
   #First mutation is updated with right values

	assert 15 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[2][0]
	assert 34 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[2][1]
	assert 49 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[2][2]


   #No other mutations, should not be updated

	assert 0 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[3][0]
	assert 0 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[3][1]
	assert 0 == DNA_bin.detecting_mutations(array.array('l',[-983172758, 17224372]))[3][2]
  

