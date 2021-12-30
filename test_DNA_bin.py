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