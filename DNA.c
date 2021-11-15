#include <Python.h>
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

	if(!PyArg_ParseTuple(args, "O", &obj))
	    return NULL;

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

	return Py_BuildValue("s", generating_mRNA(view.buf, view.shape[0]));

}

//////////////// Detecting genes
static PyObject *DNA_detecting_genes(PyObject *self, PyObject *args)
{
	Py_buffer view1;
	unsigned long long n = 0;
  	Py_buffer view3;
  	Py_buffer view4;
  	PyObject *obj1 = NULL;
  	PyObject *obj3 = NULL;
  	PyObject *obj4 = NULL;

	if(!PyArg_ParseTuple(args, "OKOO", &obj1, &n, &obj3, &obj4))
	    return NULL;

	if (PyObject_GetBuffer(obj1, &view1, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

	if (PyObject_GetBuffer(obj3, &view3, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

	if (PyObject_GetBuffer(obj4, &view4, PyBUF_ANY_CONTIGUOUS | PyBUF_FORMAT) == -1)
    	return NULL;

    if (view1.ndim != 1 || view3.ndim != 1 || view4.ndim != 1)
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array.");
		PyBuffer_Release(&view1);
		PyBuffer_Release(&view3);
		PyBuffer_Release(&view4);
		return NULL;
    }

    if (strcmp(view1.format, "I") || strcmp(view3.format, "K") || strcmp(view4.format, "K"))
    {
		PyErr_SetString(PyExc_TypeError, "Expecting a 1-dimensional array of unsigned int and unsigned long long.");
		PyBuffer_Release(&view1);
		PyBuffer_Release(&view3);
		PyBuffer_Release(&view3);
		return NULL;     
    }

	gene_map_t *gene_map;
	gene_map.genes_counter = n;
	gene_map->gene_start = view3.buf;
	gene_map->gene_end = view4;.buf

	detecting_genes(view1.buf, gene_map);
    Py_RETURN_NONE;
}

//////////////// Generating an amino acid chain (protein)
static PyObject *DNA_generating_amino_acid_chain(PyObject *self, PyObject *args)
{
	    
    Py_RETURN_NONE;
}

//////////////// Detecting probable mutation zones
static PyObject *DNA_detecting_mutations(PyObject *self, PyObject *args)
{
	Py_buffer view;
	PyObject *obj = NULL;

	if(!PyArg_ParseTuple(args, "O", &obj))
	    return NULL;

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
    return Py_BuildValue("p", detecting_mutations(view.buf, view.shape[0]));
}

//////////////// Calculating the matching score of two sequences
static PyObject *DNA_calculating_matching_score(PyObject *self, PyObject *args)
{
	    
    Py_RETURN_NONE;
}

//////////////// Hamming calculation
static PyObject *DNA_hamming(PyObject *self, PyObject *args)
{	
	int b1 = 0;
	int b2 = 0;
	if(!PyArg_ParseTuple(args, "ii", &b1, &b2))
	    return NULL;

	hamming(b1,b1);
    Py_RETURN_NONE;
}

static PyMethodDef DP_methods[] = {
	{ "generating_mRNA", DNA_generating_mRNA, METH_VARARGS, "Convert a binary DNA sequence to a string mRNA sequence"},
	{ "detecting_genes", DNA_detecting_genes, METH_VARARGS, "Detect genes"},
	{ "generating_amino_acid_chain", DNA_generating_amino_acid_chain, METH_VARARGS, "Generate an amino acid chain (protein)"},
	{ "detecting_mutations", DNA_detecting_mutations, METH_VARARGS, "Detecte probable mutation zones"},
	{ "calculating_matching_score,", DNA_calculating_matching_score, METH_VARARGS, "Detecte probable mutation zones"},
	{ "hamming", DNA_hamming, METH_VARARGS, "Hamming distance calculation"},
	{ "version", (PyCFunction)DNA_version, METH_VARARGS, "Return the version of the DNA library."},
	{NULL, NULL, 0, NULL}
};

static PyModuleDef DP_module = {
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