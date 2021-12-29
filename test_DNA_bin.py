import pytest
import DNA_bin
import array
import moduleDNA as m

###### Verify if C extension in Python is working #####

def test_get_binary_value():
	res = DNA_bin.get_binary_value(1,0)
	assert res == 1
	res = DNA_bin.get_binary_value(1,17)
	assert res == 0




