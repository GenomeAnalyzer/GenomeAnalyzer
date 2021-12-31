import pytest
import DNA_bin
import array
import moduleDNA as m

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