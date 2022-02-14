#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../headers/gene_bin.h"


/********** C-PYTHON INTERFACE DEFINTIONS **********/

static PyObject* DNAb_error = NULL;

static PyObject* DNAb_version(PyObject* self) {
	return Py_BuildValue("s", "DNA version 2.0");
}


/********** C-PYTHON INTERFACE RELATED FUNCTIONS **********/

// Returns the binary array size for long int arrays in Py_buffer.
long int DNAb_get_binary_size(Py_buffer view) {
	long int* res, res2;
	long int size = 0;
	// Copy buffer into var. memcpy needed.
	for (int i = 0; i < view.len / view.itemsize; i++) {
		res = view.buf + i * view.itemsize;
		memcpy(&res2, res, sizeof(res2));
		while (res2 > 0) {
			res2 = res2 / 2;
			size += 1;
		}
		if (size % 2 != 0) size++;
	}
	return size;
}


/********** BINARIES FUNCTION **********/

static PyObject* DNAb_get_binary_value(PyObject* self, PyObject* args) {
	Py_buffer view_seq_bin;
	PyObject* obj_seq_bin = NULL;
	int pos = 0;

	//Get the parameters (1-dimensional array of long int, int position)
	if (!PyArg_ParseTuple(args, "Oi", &obj_seq_bin, &pos))
		return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj_seq_bin, &view_seq_bin, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	if (view_seq_bin.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view_seq_bin);
		return NULL;
	}

	if (strcmp(view_seq_bin.format, "l")) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of long int.");
		PyBuffer_Release(&view_seq_bin);
		return NULL;
	}

	return Py_BuildValue("l", get_binary_value(view_seq_bin.buf, pos));
}


/******** DNA & GENES FUNCTIONS *********/

//////////////// Convert to binary
static PyObject* DNAb_convert_to_binary(PyObject* self, PyObject* args) {
	char* seq_char;
	int seq_char_size = 0;
	int seq_bin_size = 0;

	//Get the parameters (1-dimensional array of char, and its length)
	if (!PyArg_ParseTuple(args, "sii", &seq_char, &seq_char_size, &seq_bin_size))
		return NULL;

	long int* array = calloc(seq_bin_size, sizeof(long int));
	convert_to_binary(array, seq_char, seq_char_size);

	PyObject* pylist = PyList_New(seq_bin_size);

	for (int i = 0; i < seq_bin_size; ++i)
		PyList_SetItem(pylist, i, PyLong_FromLong(array[i]));

	free(array);

	return pylist;
}

//////////////// Binary to DNA
static PyObject* DNAb_binary_to_dna(PyObject* self, PyObject* args) {
	Py_buffer view_bin_dna_seq;
	PyObject* obj_bin_dna_seq = NULL;

	//Get the parameter (1-dimensional array of long int)
	if (!PyArg_ParseTuple(args, "O", &obj_bin_dna_seq))
		return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj_bin_dna_seq, &view_bin_dna_seq, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	if (view_bin_dna_seq.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view_bin_dna_seq);
		return NULL;
	}

	if (strcmp(view_bin_dna_seq.format, "l")) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of long int.");
		PyBuffer_Release(&view_bin_dna_seq);
		return NULL;
	}

	long int size = DNAb_get_binary_size(view_bin_dna_seq);

	//Return the char* value as a Python string object
	return Py_BuildValue("y", binary_to_dna(view_bin_dna_seq.buf, size));
}

//////////////// Generating mRNA
static PyObject* DNAb_generating_mRNA(PyObject* self, PyObject* args) {
	Py_buffer view_gene_seq;
	PyObject* obj_gene_seq = NULL;

	int start_pos = 0;
	int seq_size = 0;

	//Get the parameters (1-dimensional array of long int, the start position and its length)
	if (!PyArg_ParseTuple(args, "Oii", &obj_gene_seq, &start_pos, &seq_size))
		return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj_gene_seq, &view_gene_seq, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	if (view_gene_seq.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}

	if (strcmp(view_gene_seq.format, "l")) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of long int.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}

	//Return the char* value as a Python string object
	return Py_BuildValue("y", generating_mRNA(view_gene_seq.buf, start_pos, seq_size));
}

//////////////// Detecting genes
static PyObject* DNAb_detecting_genes(PyObject* self, PyObject* args) {
	Py_buffer view_gene;
	PyObject* obj_gene = NULL;

	//Get the parameter (1-dimensional arrays of long int)
	if (!PyArg_ParseTuple(args, "O", &obj_gene))
		return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj_gene, &view_gene, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	if (view_gene.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view_gene);
		return NULL;
	}

	if (strcmp(view_gene.format, "l")) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of long int");
		PyBuffer_Release(&view_gene);
		return NULL;
	}

	unsigned long gene_size = view_gene.shape[0] * 64;

	gene_map_t g;
	g.gene_start = malloc(sizeof(unsigned long long) * gene_size);
	g.gene_end = malloc(sizeof(unsigned long long) * gene_size);
	memset(g.gene_start, 0, sizeof(unsigned long long) * gene_size);
	memset(g.gene_end, 0, sizeof(unsigned long long) * gene_size);

	detecting_genes(view_gene.buf, gene_size, &g);


	PyObject* List = PyList_New(0);
	for (unsigned long long i = 0; i < g.genes_counter; i++) {
		PyObject* l = PyList_New(2);
		PyList_SET_ITEM(l, 0, PyLong_FromUnsignedLongLong(g.gene_start[i]));
		PyList_SET_ITEM(l, 1, PyLong_FromUnsignedLongLong(g.gene_end[i]));
		PyList_Append(List, l);
	}

	free(g.gene_end);
	free(g.gene_start);

	return List;
}

//////////////// Generating an amino acid chain (protein)
static PyObject* DNAb_generating_amino_acid_chain(PyObject* self, PyObject* args) {
	Py_buffer view_gene_seq;
	PyObject* obj_gene_seq = NULL;

	int start_pos = 0, seq_size = 0;

	//Get the parameters (1-dimensional array of long int, the start position, and its length)
	if (!PyArg_ParseTuple(args, "Oii", &obj_gene_seq, &start_pos, &seq_size))
		return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj_gene_seq, &view_gene_seq, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	if (view_gene_seq.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}

	if (strcmp(view_gene_seq.format, "l")) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of long int.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}

	//Return the char* value as a Python string object
	return Py_BuildValue("y", generating_amino_acid_chain(view_gene_seq.buf, start_pos, seq_size));
}

//////////////// Detecting probable mutation zones
static PyObject* DNAb_detecting_mutations(PyObject* self, PyObject* args) {
	Py_buffer view_gene_seq;
	PyObject* obj_gene_seq = NULL;

	long start_pos = 0, size_sequence = 0;

	//Get the parameters (1-dimensional arrays of long int, the start position, and its length)
	if (!PyArg_ParseTuple(args, "Oii", &obj_gene_seq, &start_pos, &size_sequence))
		return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj_gene_seq, &view_gene_seq, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	if (view_gene_seq.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}

	if (strcmp(view_gene_seq.format, "l")) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of long int.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}

	mutation_map m;

	m.size = aligned_alloc(sizeof(unsigned long), sizeof(unsigned long) * 5);
	m.start_mut = aligned_alloc(sizeof(unsigned long), sizeof(unsigned long) * 5);
	m.end_mut = aligned_alloc(sizeof(unsigned long), sizeof(unsigned long) * 5);

	//Initializing to 0 
	for (int i = 0; i < 5; i++) {
		m.size[i] = 0;
		m.start_mut[i] = 0;
		m.end_mut[i] = 0;
	}

	detecting_mutations(view_gene_seq.buf, start_pos, size_sequence, m);

	PyObject* List = PyList_New(0);
	for (short int i = 0; i < 5; i++) {
		if (m.size[i] == 0)
			continue;
		PyObject* l = PyList_New(3);
		PyList_SET_ITEM(l, 0, PyLong_FromUnsignedLong(m.size[i]));
		PyList_SET_ITEM(l, 1, PyLong_FromUnsignedLong(m.start_mut[i]));
		PyList_SET_ITEM(l, 2, PyLong_FromUnsignedLong(m.end_mut[i]));
		PyList_Append(List, l);
	}

	free(m.end_mut);
	free(m.size);
	free(m.start_mut);

	return List;
}

//////////////// Calculating the matching score of two sequences
static PyObject* DNAb_calculating_matching_score(PyObject* self, PyObject* args) {
	Py_buffer view_seq_bin1;
	Py_buffer view_seq_bin2;
	PyObject* obj_seq_bin1 = NULL;
	PyObject* obj_seq_bin2 = NULL;

	long start_pos1 = 0, start_pos2 = 0, seq_size1 = 0, seq_size2 = 0;

	//Get the parameters (2 1-dimensional arrays, with its start position and its length)
	if (!PyArg_ParseTuple(args, "OiiOii", &obj_seq_bin1, &start_pos1, &seq_size1, &obj_seq_bin2, &start_pos2, &seq_size2))
		return NULL;

	//Get the first array memory view
	if (PyObject_GetBuffer(obj_seq_bin1, &view_seq_bin1, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	//Get the second array memory view
	if (PyObject_GetBuffer(obj_seq_bin2, &view_seq_bin2, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	if (view_seq_bin1.ndim != 1 || view_seq_bin2.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "Expecting 2 1-dimensional array.");
		PyBuffer_Release(&view_seq_bin1);
		PyBuffer_Release(&view_seq_bin2);
		return NULL;
	}

	if (strcmp(view_seq_bin1.format, "l") || strcmp(view_seq_bin2.format, "l")) {
		PyErr_SetString(PyExc_TypeError, "Expecting 2 1-dimensional array of long int.");
		PyBuffer_Release(&view_seq_bin1);
		PyBuffer_Release(&view_seq_bin2);
		return NULL;
	}

	//Return the float value as a Python float object
	return Py_BuildValue("f", calculating_matching_score(view_seq_bin1.buf, start_pos1, seq_size1, view_seq_bin2.buf, start_pos2, seq_size2));
}


/********** C-PYTHON INTERFACE SETUP FUNCTIONS **********/

//Register the methods to be made available Python side
static PyMethodDef DNAb_methods [] = {
	{ "get_binary_value", DNAb_get_binary_value, METH_VARARGS, "Retrieve one bit from the binary array sequence"},
	{ "convert_to_binary", DNAb_convert_to_binary, METH_VARARGS, "Convert a DNA base sequence to its binary array format"},
	{ "binary_to_dna", DNAb_binary_to_dna, METH_VARARGS, "Convert a DNA sequence in binary array format to its DNA bases"},
	{ "generating_mRNA", DNAb_generating_mRNA, METH_VARARGS, "Convert a DNA sequence in binary array format to its mRNA sequence"},
	{ "detecting_genes", DNAb_detecting_genes, METH_VARARGS, "Detects genes in the mRNA sequence in binary array format and maps them"},
	{ "generating_amino_acid_chain", DNAb_generating_amino_acid_chain, METH_VARARGS, "Generate an amino acid chain (protein) from a binary arary sequence"},
	{ "detecting_mutations", DNAb_detecting_mutations, METH_VARARGS, "Detects probable mutation areas"},
	{ "calculating_matching_score", DNAb_calculating_matching_score, METH_VARARGS, "Calculates the matching score of two binary array sequences"},
	{ "version", (PyCFunction)DNAb_version, METH_VARARGS, "Return the version of the DNA library"},
	{NULL, NULL, 0, NULL}
};

// Define the python module
static PyModuleDef DNAb_module = {
	PyModuleDef_HEAD_INIT,
	"DNA_bin",
	"Genome analyser",
	-1,
	DNAb_methods
};

PyMODINIT_FUNC PyInit_DNA_bin() {
	PyObject* obj = PyModule_Create(&DNAb_module);

	if (!obj)
		return NULL;

	DNAb_error = PyErr_NewException("DNAb.error", NULL, NULL);
	Py_XINCREF(DNAb_error);

	if (PyModule_AddObject(obj, "error", DNAb_error) < 0) {
		Py_XDECREF(DNAb_error);
		Py_CLEAR(DNAb_error);
		Py_DECREF(obj);
		return NULL;
	}

	return obj;
}
