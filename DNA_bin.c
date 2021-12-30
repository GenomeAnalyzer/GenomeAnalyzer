#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gene_bin.h"


static PyObject *DNAb_error = NULL;

static PyObject *DNAb_version(PyObject *self)
{
	return Py_BuildValue("s", "DNA version 2.0");
}

/********** BINARIES FUNCTION **********/

static PyObject *DNAb_get_binary_value(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;
	int obj2 = 0;

	//Get the parameter (1-dimensional array of unsigned int, int position)
	if(!PyArg_ParseTuple(args, "Oi", &obj, &obj2))
	    return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj, &view, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

    if (view.ndim != 1)
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view);
		return NULL;
    }

    if (strcmp(view.format, "I"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned int.");
		PyBuffer_Release(&view);
		return NULL;     
    }

	return Py_BuildValue("i", get_binary_value(view.buf, obj2));
}

static PyObject *DNAb_change_binary_value(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;
	int obj2 = 0;
	int obj3 = 0;

	//Get the parameter (1-dimensional array of unsigned int, int position, int value)
	if(!PyArg_ParseTuple(args, "Oii", &obj, &obj2, &obj3))
	    return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj, &view, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

    if (view.ndim != 1)
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view);
		return NULL;
    }

    if (strcmp(view.format, "I"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned int.");
		PyBuffer_Release(&view);
		return NULL;     
    }

	unsigned int *array = change_binary_value(view.buf, obj2, obj3);
	PyObject *pylist = PyList_New(view.shape[0]);

	for (int i = 0; i < view.shape[0]; i++)
		PyList_SetItem(pylist, i, PyLong_FromLong(array[i]));

	return pylist;
}

static PyObject *DNAb_set_binary_array(PyObject *self, PyObject *args)
{
	// unsigned int* set_binary_array(const char *seq_char, const unsigned seq_size){
	Py_buffer view;
	PyObject* obj = NULL;
	char *seq_char;
	int seq_size;

	//Get the parameter (1-dimensional array of char, and its length)
	// if (!PyArg_ParseTuple(args, "y#", &seq_char, &seq_size))
	if (!PyArg_ParseTuple(args, "s#", &seq_char, &seq_size))
		return NULL;

	printf("ok here\n");
	printf("seq_size : %d\n", seq_size);
	printf("seq_char : %p\n", (void*) seq_char);
	// TODO : corriger le segfault sur seq_char
	// printf("seq_char : %s\n", *seq_char);
	// printf("seq_char size: %d\n", sizeof(seq_char));

	unsigned int* array = set_binary_array(seq_char, seq_size);
	PyObject* pylist = PyList_New(seq_size*2);

	for (int i = 0; i < seq_size*2; i++)
		PyList_SetItem(pylist, i, PyLong_FromLong(array[i]));

	return pylist;
}

static PyObject *DNAb_xor_binary_array(PyObject *self, PyObject *args)
{
	
}

static PyObject *DNAb_popcount_binary_array(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;

	//Get the parameter (1-dimensional array of unsigned int)
	if(!PyArg_ParseTuple(args, "O", &obj))
	    return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj, &view, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

    if (view.ndim != 1)
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view);
		return NULL;
    }

    if (strcmp(view.format, "I"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned int.");
		PyBuffer_Release(&view);
		return NULL;     
    }

	// Copy buffer into var. memcpy needed otherwise it doesn't work.
	unsigned int *res = view.buf;
	unsigned int res2;
	memcpy(&res2, res, sizeof(res));
	unsigned int size = 0;
	// Calculate binary size.
	while (res2 > 0){
		res2 = res2 / 2;
		size += 1;
	}

	return Py_BuildValue("i", popcount_binary_array(view.buf, size));
}



/******** DNA & GENES FUNCTION *********/

//////////////// Convert to binary
static PyObject *DNAb_convert_to_binary(PyObject *self, PyObject *args)
{
	char *obj = NULL;
	unsigned obj2 = 0;

	//Get the parameter (char* value)
	if(!PyArg_ParseTuple(args, "si", &obj,&obj2))
	    return NULL;

	obj2 *=2;
	unsigned int *test = convert_to_binary(obj,obj2);
	PyObject *pylist = PyList_New(obj2);

	for (unsigned  i=0; i < obj2; i++)
		PyList_SetItem(pylist, i, PyLong_FromLong(test[i]));

	return pylist;
}

//////////////// Binary to DNA
static PyObject *DNAb_binary_to_dna(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;

	//Get the parameter (1-dimensional array of unsigned int)
	if(!PyArg_ParseTuple(args, "O", &obj))
	    return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj, &view, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

    if (view.ndim != 1)
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view);
		return NULL;
    }

    if (strcmp(view.format, "I"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned int.");
		PyBuffer_Release(&view);
		return NULL;     
    }

    //Return the char* value as a Python string object
	return Py_BuildValue("y", generating_mRNA(view.buf, view.shape[0]));
}

//////////////// Generating mRNA
static PyObject *DNAb_generating_mRNA(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;

	//Get the parameter (1-dimensional array of unsigned int)
	if(!PyArg_ParseTuple(args, "O", &obj))
	    return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj, &view, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

    if (view.ndim != 1)
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view);
		return NULL;
    }

    if (strcmp(view.format, "I"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned int.");
		PyBuffer_Release(&view);
		return NULL;     
    }

    //Return the char* value as a Python string object
	return Py_BuildValue("y", generating_mRNA(view.buf, view.shape[0]));
}

//////////////// Detecting genes
static PyObject *DNAb_detecting_genes(PyObject *self, PyObject *args)
{
	Py_buffer view;
  	PyObject *obj = NULL;

  	//Get the parameter (1-dimensional arrays of unsigned int)
	if (!PyArg_ParseTuple(args, "O", &obj))
	    return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj, &view, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

    if (view.ndim != 1)
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view);
		return NULL;
    }

    if (strcmp(view.format, "I"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned int");
		PyBuffer_Release(&view);
		return NULL;     
    }

	gene_map_t g;

    g.gene_start = malloc(sizeof(*g.gene_start) * view.shape[0]);
    g.gene_end = malloc(sizeof(*g.gene_end) * view.shape[0]);

	detecting_genes(view.buf, view.shape[0], &g);

	PyObject *List = PyList_New(0);
	for(unsigned long long i = 0; i < g.genes_counter; i ++)
    {
		PyObject *l = PyList_New(2);
    	PyList_SET_ITEM(l, 0, PyLong_FromUnsignedLongLong(g.gene_start[i]));
    	PyList_SET_ITEM(l, 1, PyLong_FromUnsignedLongLong(g.gene_end[i]));
    	PyList_Append(List, l);    	
    }

	return List;
}

//////////////// Generating an amino acid chain (protein)
static PyObject *DNAb_generating_amino_acid_chain(PyObject *self, PyObject *args)
{
    Py_buffer view;
	PyObject *obj = NULL;

	//Get the parameter (1-dimensional array of unsigned int)
	if(!PyArg_ParseTuple(args, "O", &obj))
	    return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj, &view, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
        return NULL;

    if (view.ndim != 1)
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view);
		return NULL;
    }

    if (strcmp(view.format, "I"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned int.");
		PyBuffer_Release(&view);
		return NULL;
    }

    //Return the char* value as a Python string object
	return Py_BuildValue("y", generating_amino_acid_chain(view.buf, view.shape[0]));
}

//////////////// Detecting probable mutation zones
static PyObject *DNAb_detecting_mutations(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;

	//Get the parameter (1-dimensional arrays of unsigned int)
	if(!PyArg_ParseTuple(args, "O", &obj))
	    return NULL;

	//Get the array memory view
	if (PyObject_GetBuffer(obj, &view, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

    if (view.ndim != 1)
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view);
		return NULL;
    }

    if (strcmp(view.format, "I"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned int.");
		PyBuffer_Release(&view);
		return NULL;     
    }

	mutation_map m;

	m.size = malloc(sizeof(unsigned long) * 5);
	m.start_mut = malloc(sizeof(unsigned long) * 5);
	m.end_mut = malloc(sizeof(unsigned long) * 5);

	//Initializing to 0 
	for(int i = 0; i < 5; i ++){
		m.size[i]=0;
    	m.start_mut[i]=0;
    	m.end_mut[i]=0;   	
    }

    detecting_mutations(view.buf, view.shape[0], m);

    PyObject* List = PyList_New(0);
    for(short int i = 0; i < 5; i ++){
		PyObject *l = PyList_New(3);
		PyList_SET_ITEM(l, 0, PyLong_FromUnsignedLong(m.size[i]));
    	PyList_SET_ITEM(l, 1, PyLong_FromUnsignedLong(m.start_mut[i]));
    	PyList_SET_ITEM(l, 2, PyLong_FromUnsignedLong(m.end_mut[i]));
    	PyList_Append(List, l);
    }

	return List;
}

//////////////// Calculating the matching score of two sequences
static PyObject *DNAb_calculating_matching_score(PyObject *self, PyObject *args)
{
	Py_buffer view1;
	Py_buffer view2;
	PyObject *obj1 = NULL;
	PyObject *obj2 = NULL;

	//Get the parameter (2 1-dimensional arrays)
	if(!PyArg_ParseTuple(args, "OO", &obj1, &obj2))
	    return NULL;

	//Get the first array memory view
	if (PyObject_GetBuffer(obj1, &view1, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

    //Get the second array memory view
	if (PyObject_GetBuffer(obj2, &view2, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

    if (view1.ndim != 1 || view2.ndim != 1){
		PyErr_SetString(PyExc_TypeError, "Expecting 2 1-dimensional array.");
		PyBuffer_Release(&view1);
		PyBuffer_Release(&view2);
		return NULL;
    }

    if (strcmp(view1.format, "I") || strcmp(view2.format, "I")){
		PyErr_SetString(PyExc_TypeError, "Expecting 2 1-dimensional array of unsigned int.");
		PyBuffer_Release(&view1);
		PyBuffer_Release(&view2);
		return NULL;
    }

	// printf("%ld\n", view1.shape[0]);

	//Return the float value as a Python float object
	return Py_BuildValue("f", calculating_matching_score(view1.buf, view1.shape[0], view2.buf, view2.shape[0]));
}

//Register the methods to be made available Python side
static PyMethodDef DNAb_methods[] = {
	{ "get_binary_value", DNAb_get_binary_value, METH_VARARGS, ""},	
	{ "change_binary_value", DNAb_change_binary_value, METH_VARARGS, ""},	
	{ "set_binary_array", DNAb_set_binary_array, METH_VARARGS, ""},	
	{ "xor_binary_array", DNAb_xor_binary_array, METH_VARARGS, ""},	
	{ "popcount_binary_array", DNAb_popcount_binary_array, METH_VARARGS, ""},	
	{ "convert_to_binary", DNAb_convert_to_binary, METH_VARARGS, "Convert a char sequence to a binary sequence"},	
	{ "binary_to_dna", DNAb_binary_to_dna, METH_VARARGS, "Convert a binary sequence to a DNA sequence"},	
	{ "generating_mRNA", DNAb_generating_mRNA, METH_VARARGS, "Convert a binary DNA sequence to a string mRNA sequence"},
	{ "detecting_genes", DNAb_detecting_genes, METH_VARARGS, "Detect genes"},
	{ "generating_amino_acid_chain", DNAb_generating_amino_acid_chain, METH_VARARGS, "Generate an amino acid chain (protein)"},
	{ "detecting_mutations", DNAb_detecting_mutations, METH_VARARGS, "Detecte probable mutation zones"},
	{ "calculating_matching_score", DNAb_calculating_matching_score, METH_VARARGS, "Detecte probable mutation zones"},
	{ "version", (PyCFunction)DNAb_version, METH_VARARGS, "Return the version of the DNA library."},
	{NULL, NULL, 0, NULL}
};

static PyModuleDef DNAb_module = {
	PyModuleDef_HEAD_INIT,
	"DNA_bin",
	"Genome analyser",
	-1,
	DNAb_methods
};

PyMODINIT_FUNC PyInit_DNA_bin()
{
	PyObject *obj = PyModule_Create(&DNAb_module);

  	if (!obj)
    	return NULL;

  	DNAb_error = PyErr_NewException("DNAb.error", NULL, NULL);
  	Py_XINCREF(DNAb_error);
  
  	if (PyModule_AddObject(obj, "error", DNAb_error) < 0)
    {
      	Py_XDECREF(DNAb_error);
      	Py_CLEAR(DNAb_error);
      	Py_DECREF(obj);
      	return NULL;
    }

  return obj;
}