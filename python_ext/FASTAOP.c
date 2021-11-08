#include <Python.h>
#include "FASTAOP_core.h"

static PyObject* FASTAOP_error = NULL;

// http://web.mit.edu/people/amliu/vrut/python/ext/


/// Functions

static PyObject* FASTAOP_version(PyObject* self) {
    return Py_BuildValue("s", "FASTAOP version 0.1");
}

static PyObject* FASTAOP_hamming(PyObject* self, PyObject* args) {
    int b1 = 0;
    int b2 = 0;
    if (!PyArg_ParseTuple(args, "ii", &b1, &b2)) return NULL;

    return Py_BuildValue("i", hamming(b1, b2));
}

static PyObject* FASTAOP_mRNA_generate(PyObject* self, PyObject* args) {
    char* s;
    int size;
    if (!PyArg_ParseTuple(args, "s#", &s, &size)) return NULL;

    return Py_BuildValue("s", RNA_generate(size, s));
}

static PyObject* FASTAOP_genes_detect(PyObject* self, PyObject* args) {
    int var;
    if (!PyArg_ParseTuple(args, "i", &var)) return NULL;
    return Py_BuildValue("i", genes_detect(var));
}

static PyObject* FASTAOP_amino_acid_chain_generate(PyObject* self, PyObject* args) {
    int var;
    if (!PyArg_ParseTuple(args, "i", &var)) return NULL;
    return Py_BuildValue("i", amino_acid_chain_generate(var));
}

static PyObject* FASTAOP_mutations_detect(PyObject* self, PyObject* args) {
    int var;
    if (!PyArg_ParseTuple(args, "i", &var)) return NULL;
    return Py_BuildValue("i", mutations_detect(var));
}

static PyObject* FASTAOP_matching_score(PyObject* self, PyObject* args) {
    int a, b;
    if (!PyArg_ParseTuple(args, "ii", &a, &b)) return NULL;

    return Py_BuildValue("i", matching_score_calculate(a, b));
}


/// Register the methods to be made available Python side
static PyMethodDef FASTAOP_methods [] = {
    {"hamming", FASTAOP_hamming, METH_VARARGS, "Calculation Hamming distance"},
    {"RNA_generate", FASTAOP_mRNA_generate, METH_VARARGS, "Generating mRNA"},
    {"genes_detect", FASTAOP_genes_detect, METH_VARARGS, "Detecting genes"},
    {"amino_acid_chain_generate", FASTAOP_amino_acid_chain_generate, METH_VARARGS, "Generating amino acid chain"},
    {"mutations_detect", FASTAOP_mutations_detect, METH_VARARGS, "Detecting mutation"},
    {"matching_score_calculate", FASTAOP_matching_score, METH_VARARGS, "Calculating matching score"},

    {"version", (PyCFunction)FASTAOP_version, METH_VARARGS, "Returns the version of the FASTAOP library."},

    {NULL, NULL, 0, NULL}
};

static PyModuleDef FASTAOP_module = {
    PyModuleDef_HEAD_INIT,
    "FASTAOP",
    "FASTAOP Basic FASTA Operations Library",
    -1,
    FASTAOP_methods
};

///
PyMODINIT_FUNC PyInit_FASTAOP() {
    PyObject* obj = PyModule_Create(&FASTAOP_module);

    if (!obj)
        return NULL;

    FASTAOP_error = PyErr_NewException("FASTAOP.error", NULL, NULL);
    Py_XINCREF(FASTAOP_error);

    if (PyModule_AddObject(obj, "error", FASTAOP_error) < 0)     {
        Py_XDECREF(FASTAOP_error);
        Py_CLEAR(FASTAOP_error);
        Py_DECREF(obj);
        return NULL;
    }

    return obj;
}
