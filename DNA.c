#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gene.h"

static PyObject *DNA_error = NULL;

static PyObject *DNA_version(PyObject *self)
{
	return Py_BuildValue("s", "DNA version 0.1");
}

//////////////// Convert to binary
static PyObject *DNA_convert_to_binary(PyObject *self, PyObject *args)
{
	char *obj = NULL;

	unsigned obj2 = 0;
	//Get the parameter (char* value)
	if(!PyArg_ParseTuple(args, "si", &obj,&obj2))
	    return NULL;

	obj2 *=2;
	unsigned short *test = convert_to_binary(obj,obj2);
	 PyObject *pylist, *item;
	pylist = PyList_New(obj2);

	for (unsigned  i=0; i < obj2; i++) {
		item = PyLong_FromLong(test[i]);
		PyList_SetItem(pylist, i, item);
	}
return pylist;
}

//////////////// Generating mRNA
static PyObject *DNA_generating_mRNA(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;

	//Get the parameter (1-dimensional array of unsigned short)
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

    if (strcmp(view.format, "H"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned short.");
		PyBuffer_Release(&view);
		return NULL;     
    }

    //Return the char* value as a Python string object
	return Py_BuildValue("y", generating_mRNA(view.buf, view.shape[0]));
}

//////////////// Detecting genes
static PyObject *DNA_detecting_genes(PyObject *self, PyObject *args)
{
	Py_buffer view;
  	PyObject *obj = NULL;

  	//Get the parameter (1-dimensional arrays of unsigned short)
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

    if (strcmp(view.format, "L"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned long");
		PyBuffer_Release(&view);
		return NULL;     
    }

	gene_map_t g;

    g.gene_start = malloc(sizeof(*g.gene_start) * view.shape[0]+1);
    g.gene_end = malloc(sizeof(*g.gene_end) * view.shape[0]+1);


	detecting_genes(view.buf, view.shape[0], &g);


	PyObject *List = PyList_New(0);
	for(unsigned long long i = 0; i < g.genes_counter; i ++)
    {
		PyObject *l = PyList_New(2);
		//printf("%llu truc la \n",g.gene_start[i]);
    	PyList_SET_ITEM(l, 0, PyLong_FromUnsignedLongLong(g.gene_start[i]));
    	PyList_SET_ITEM(l, 1, PyLong_FromUnsignedLongLong(g.gene_end[i]));
    	PyList_Append(List, l);    	
    }

	free(g.gene_end);
	free(g.gene_start);

	return List;
}

//////////////// Generating an amino acid chain (protein)
static PyObject *DNA_generating_amino_acid_chain(PyObject *self, PyObject *args)
{
    Py_buffer view;
	PyObject *obj = NULL;

	//Get the parameter (1-dimensional array of unsigned short)
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

    if (strcmp(view.format, "H"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned short.");
		PyBuffer_Release(&view);
		return NULL;
    }

    //Return the char* value as a Python string object
	return Py_BuildValue("y", generating_amino_acid_chain(view.buf, view.shape[0]));
}

//////////////// Detecting probable mutation zones
static PyObject *DNA_detecting_mutations(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;

	//Get the parameter (1-dimensional arrays of unsigned short)
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

    if (strcmp(view.format, "H"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned short.");
		PyBuffer_Release(&view);
		return NULL;     
    }

	mutation_map m;

	m.size = aligned_alloc(sizeof(unsigned long), sizeof(unsigned long)*5);
	m.start_mut = aligned_alloc(sizeof(unsigned long), sizeof(unsigned long)*5);
	m.end_mut = aligned_alloc(sizeof(unsigned long), sizeof(unsigned long)*5);
	//Initializing to 0 
	for(int i = 0; i < 5; i++){
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

    free(m.size);
    free(m.start_mut);
    free(m.end_mut);

	return List;
}

//////////////// Calculating the matching score of two sequences
static PyObject *DNA_calculating_matching_score(PyObject *self, PyObject *args)
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

    if (strcmp(view1.format, "H") || strcmp(view2.format, "H")){
		PyErr_SetString(PyExc_TypeError, "Expecting 2 1-dimensional array of unsigned short.");
		PyBuffer_Release(&view1);
		PyBuffer_Release(&view2);
		return NULL;
    }

	// printf("%ld\n", view1.shape[0]);

	//Return the float value as a Python float object
	return Py_BuildValue("f", calculating_matching_score(view1.buf, view1.shape[0], view2.buf, view2.shape[0]));
}

//Register the methods to be made available Python side
static PyMethodDef DNA_methods[] = {
	{ "convert_to_binary", DNA_convert_to_binary, METH_VARARGS, "Convert a char sequence to a binary sequence"},	
	{ "generating_mRNA", DNA_generating_mRNA, METH_VARARGS, "Convert a binary DNA sequence to a string mRNA sequence"},
	{ "detecting_genes", DNA_detecting_genes, METH_VARARGS, "Detect genes"},
	{ "generating_amino_acid_chain", DNA_generating_amino_acid_chain, METH_VARARGS, "Generate an amino acid chain (protein)"},
	{ "detecting_mutations", DNA_detecting_mutations, METH_VARARGS, "Detecte probable mutation zones"},
	{ "calculating_matching_score", DNA_calculating_matching_score, METH_VARARGS, "Detecte probable mutation zones"},
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