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