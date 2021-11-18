import pytest
import DNA
import array
import numpy as np

###### Verify if C extension in Python is working #####

def test_DNA_generating_mRNA():
	res = DNA.generating_mRNA(array.array('i', [0,1,0,0,1,1,0,0,0,1,1,0]))
	assert res == "GAUAGC"

# def test_DNA_detecting_genes():
# 	DNA.detecting_genes(array.array('I', [0,1,0,0,1,1,0,0,0,1,1,0]))
# 	assert True

def test_DNA_generating_amino_acid_chain():
	res = DNA.generating_amino_acid_chain("GAUAGC")
	assert res == "GAUAGC"

def test_DNA_detecting_mutations():
	res = DNA.detecting_mutations(array.array('I', [0,1,0,0,1,1,0,0,0,1,1,0]))
	assert res

# def test_DNA_calculating_matching_score():
# 	res = DNA.calculating_matching_score(0b010011000110, 0b101100111001)
# 	assert res == -1.0

def test_DNA_hamming():
	res = DNA.hamming(0b010011000110, 0b101100111001)
	assert res == 12

