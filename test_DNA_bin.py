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
	assert 0 == 0
	# Test if the algorithm is OK

	seq_char = "ACGT" # for s# objects
	seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGA"#.encode('utf-8') # for y# objects
	seq_size = 4
	
	resbin = DNA_bin.set_binary_array(seq_char)
	# TODO : corriger le segfault sur seq_char
	assert resbin == [2101911378,172292753,-265825143] 

	#char *seq_char = "GACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGAGACCTTCGA";
	#unsigned seq_size = 45;

	# unsigned int *seq_bin = NULL;
	# seq_bin = set_binary_array(seq_char, seq_size);

	# unsigned int seq_sol[3] = {2101911378, 172292753, -265825143};

	# for (int i = 0; i < 3; ++i)
	# 	assert_int_equal(seq_bin[i], seq_sol[i]);
	 
	# free(seq_bin);


#def test_popcount_binary_array():
	# Test for all binaries from 00000 to 11111
#	for i in range(0,32):
#		popc_expected_result = int(i%2) + int(i/2%2) + int(i/4%2) + int(i/8%2) + int(i/16%2)
#		popc_result = DNA_bin.popcount_binary_array(array.array('I', [i]))
#		assert popc_expected_result == popc_result

	# 9350764 = 001101100111010101110001
	#assert 13 == DNA_bin.popcount_binary_array(array.array('I', [9350764]))