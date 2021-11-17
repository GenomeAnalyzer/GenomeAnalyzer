#include <Python.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gene.h"

static PyObject *DNA_error = NULL;

static PyObject *DNA_version(PyObject *self)
{
	return Py_BuildValue("s", "DNA version 0.1");
}

//////////////// Generating mRNA
static PyObject *DNA_generating_mRNA(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;

	//Get the parameter (1-dimensional arrays)
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

    if (strcmp(view.format, "i"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of floats.");
		PyBuffer_Release(&view);
		return NULL;     
    }

    //Return the char* value as a Python string object
	return Py_BuildValue("s", generating_mRNA(view.buf, view.shape[0]));

}

//////////////// Detecting genes
static PyObject *DNA_detecting_genes(PyObject *self, PyObject *args)
{
	Py_buffer view;
  	PyObject *obj = NULL;

  	//Get the parameter (1-dimensional arrays)
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
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned int");
		PyBuffer_Release(&view);
		return NULL;     
    }

	gene_map_t g;
	g.genes_counter = 0;
	g.gene_start = malloc(sizeof(int*) * MAX_GENES);
	g.gene_end = malloc(sizeof(int*) * MAX_GENES);

	detecting_genes(view.buf, &g);

	//Return nothing because void type
    Py_RETURN_NONE;
}

//////////////// Generating an amino acid chain (protein)
static PyObject *DNA_generating_amino_acid_chain(PyObject *self, PyObject *args)
{
	char *obj = NULL;

	//Get the parameter (char* value)
	if(!PyArg_ParseTuple(args, "s", &obj))
	    return NULL;

	//Return the char* value as a Python string object
    return Py_BuildValue("s", generating_amino_acid_chain(obj));
}

//////////////// Detecting probable mutation zones
static PyObject *DNA_detecting_mutations(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;

	//Get the parameter (1-dimensional arrays)
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

    //Return the boolean value as a Python boolean object
    return Py_BuildValue("O", detecting_mutations(view.buf, view.shape[0]) ? Py_True : Py_False);
}

//////////////// Calculating the matching score of two sequences
static PyObject *DNA_calculating_matching_score(PyObject *self, PyObject *args)
{
	int obj1 = 0;
	int obj2 = 0;

	//Get the parameter (2 int values)
	if(!PyArg_ParseTuple(args, "ii", &obj1, &obj2))
		return NULL;

	//Return the float value as a Python float object
    return Py_BuildValue("f", calculating_matching_score(obj1, obj2));
}

//////////////// Hamming calculation
static PyObject *DNA_hamming(PyObject *self, PyObject *args)
{	
	int obj1 = 0;
	int obj2 = 0;

	//Get the parameter (2 int values)
	if(!PyArg_ParseTuple(args, "ii", &obj1, &obj2))
	    return NULL;

	//Return the int value as a Python int object
    return Py_BuildValue("i", hamming(obj1, obj2));
}

//Register the methods to be made available Python side
static PyMethodDef DNA_methods[] = {
	{ "generating_mRNA", DNA_generating_mRNA, METH_VARARGS, "Convert a binary DNA sequence to a string mRNA sequence"},
	{ "detecting_genes", DNA_detecting_genes, METH_VARARGS, "Detect genes"},
	{ "generating_amino_acid_chain", DNA_generating_amino_acid_chain, METH_VARARGS, "Generate an amino acid chain (protein)"},
	{ "detecting_mutations", DNA_detecting_mutations, METH_VARARGS, "Detecte probable mutation zones"},
	{ "calculating_matching_score", DNA_calculating_matching_score, METH_VARARGS, "Detecte probable mutation zones"},
	{ "hamming", DNA_hamming, METH_VARARGS, "Hamming distance calculation"},
	{ "version", (PyCFunction)DNA_version, METH_VARARGS, "Return the version of the DNA library."},
	{NULL, NULL, 0, NULL}
};

static PyModuleDef DNA_module = {
	PyModuleDef_HEAD_INIT,
	"DNA",
	"Genome analyser",
	-1,
	DNA_methods
};

PyMODINIT_FUNC PyInit_DNA()
{
	PyObject *obj = PyModule_Create(&DNA_module);

  	if (!obj)
    	return NULL;

  	DNA_error = PyErr_NewException("DNA.error", NULL, NULL);
  	Py_XINCREF(DNA_error);
  
  	if (PyModule_AddObject(obj, "error", DNA_error) < 0)
    {
      	Py_XDECREF(DNA_error);
      	Py_CLEAR(DNA_error);
      	Py_DECREF(obj);
      	return NULL;
    }

  return obj;
}