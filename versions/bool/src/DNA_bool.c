#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../headers/gene_bool.h"


/********** C-PYTHON INTERFACE DEFINTIONS **********/

static PyObject* DNA_bool_error = NULL;

static PyObject* DNA_bool_version(PyObject* self) {
	return Py_BuildValue("s", "DNA version 2.0");
}




/******** DNA & GENES FUNCTIONS *********/

//////////////// Convert to binary
static PyObject* DNAb_convert_to_binary(PyObject* self, PyObject* args) {
	char* seq_char;
	long int seq_size = 0;

	//Get the parameters (1-dimensional array of char, and its length)
	if (!PyArg_ParseTuple(args, "si", &seq_char, &seq_size))
		return NULL;

	_Bool* array = convert_to_binary(seq_char, seq_size);


	PyObject* pylist = PyList_New(seq_size*2);

	for (long int i = 0; i < seq_size*2; i++)
		PyList_SetItem(pylist, i, PyBool_FromLong(array[i]));

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

	if (strcmp(view_bin_dna_seq.format, "B")) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of long int.");
		PyBuffer_Release(&view_bin_dna_seq);
		return NULL;
	}


	//Return the char* value as a Python string object
	return Py_BuildValue("y", binary_to_dna(view_bin_dna_seq.buf, view_bin_dna_seq.shape[0]));
}

//////////////// Generating mRNA
static PyObject* DNAb_generating_mRNA(PyObject* self, PyObject* args) {
	Py_buffer view_gene_seq;
	PyObject* obj_gene_seq = NULL;


	//Get the parameters (1-dimensional array of long int, the start position and its length)
	if (!PyArg_ParseTuple(args, "O", &obj_gene_seq))
		return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj_gene_seq, &view_gene_seq, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	if (view_gene_seq.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}



	if (strcmp(view_gene_seq.format, "B")) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of long int.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}

	//Return the char* value as a Python string object
	return Py_BuildValue("y", generating_mRNA(view_gene_seq.buf, view_gene_seq.shape[0]));
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

	if (strcmp(view_gene.format, "B")) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of long int");
		PyBuffer_Release(&view_gene);
		return NULL;
	}


	gene_map_t g;
	g.gene_start = malloc(sizeof(*g.gene_start) * view_gene.shape[0]/2);
	g.gene_end = malloc(sizeof(*g.gene_end) * view_gene.shape[0]/2);

	detecting_genes(view_gene.buf, view_gene.shape[0], &g);


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


	//Get the parameters (1-dimensional array of long int, the start position, and its length)
	if (!PyArg_ParseTuple(args, "O", &obj_gene_seq))
		return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj_gene_seq, &view_gene_seq, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	if (view_gene_seq.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}

	if (strcmp(view_gene_seq.format, "B")) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of long int.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}

	//Return the char* value as a Python string object
	return Py_BuildValue("y", generating_amino_acid_chain(view_gene_seq.buf, view_gene_seq.shape[0]));
}

//////////////// Detecting probable mutation zones
static PyObject* DNAb_detecting_mutations(PyObject* self, PyObject* args) {
	Py_buffer view_gene_seq;
	PyObject* obj_gene_seq = NULL;


	//Get the parameters (1-dimensional arrays of long int, the start position, and its length)
	if (!PyArg_ParseTuple(args, "O", &obj_gene_seq))
		return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj_gene_seq, &view_gene_seq, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
		return NULL;

	if (view_gene_seq.ndim != 1) {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view_gene_seq);
		return NULL;
	}

	if (strcmp(view_gene_seq.format, "B")) {
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

	detecting_mutations(view_gene_seq.buf, view_gene_seq.shape[0], m);

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


	//Get the parameters (2 1-dimensional arrays, with its start position and its length)
	if (!PyArg_ParseTuple(args, "OO", &obj_seq_bin1,&obj_seq_bin2))
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

	if (strcmp(view_seq_bin1.format, "B") || strcmp(view_seq_bin2.format, "B")) {
		PyErr_SetString(PyExc_TypeError, "Expecting 2 1-dimensional array of long int.");
		PyBuffer_Release(&view_seq_bin1);
		PyBuffer_Release(&view_seq_bin2);
		return NULL;
	}

	//Return the float value as a Python float object
	return Py_BuildValue("f", calculating_matching_score(view_seq_bin1.buf, view_seq_bin1.shape[0],view_seq_bin2.buf, view_seq_bin2.shape[0]));
}


/********** C-PYTHON INTERFACE SETUP FUNCTIONS **********/

//Register the methods to be made available Python side
static PyMethodDef DNA_bool_methods [] = {
	{ "convert_to_binary", DNAb_convert_to_binary, METH_VARARGS, "Convert a DNA base sequence to its binary array format"},
	{ "binary_to_dna", DNAb_binary_to_dna, METH_VARARGS, "Convert a DNA sequence in binary array format to its DNA bases"},
	{ "generating_mRNA", DNAb_generating_mRNA, METH_VARARGS, "Convert a DNA sequence in binary array format to its mRNA sequence"},
	{ "detecting_genes", DNAb_detecting_genes, METH_VARARGS, "Detects genes in the mRNA sequence in binary array format and maps them"},
	{ "generating_amino_acid_chain", DNAb_generating_amino_acid_chain, METH_VARARGS, "Generate an amino acid chain (protein) from a binary arary sequence"},
	{ "detecting_mutations", DNAb_detecting_mutations, METH_VARARGS, "Detects probable mutation areas"},
	{ "calculating_matching_score", DNAb_calculating_matching_score, METH_VARARGS, "Calculates the matching score of two binary array sequences"},
	{ "version", (PyCFunction)DNA_bool_version, METH_VARARGS, "Return the version of the DNA library"},
	{NULL, NULL, 0, NULL}
};

// Define the python module
static PyModuleDef DNA_bool_module = {
	PyModuleDef_HEAD_INIT,
	"DNA_boolt",
	"Genome analyser",
	-1,
	DNA_bool_methods
};

PyMODINIT_FUNC PyInit_DNA_bool() {
	PyObject* obj = PyModule_Create(&DNA_bool_module);

	if (!obj)
		return NULL;

	DNA_bool_error = PyErr_NewException("DNA_bool.error", NULL, NULL);
	Py_XINCREF(DNA_bool_error);

	if (PyModule_AddObject(obj, "error", DNA_bool_error) < 0) {
		Py_XDECREF(DNA_bool_error);
		Py_CLEAR(DNA_bool_error);
		Py_DECREF(obj);
		return NULL;
	}

	return obj;
}
