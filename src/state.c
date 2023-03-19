#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

static PyObject *view(PyObject *self, PyObject *args, PyObject *kwds)
{
  int N,K;
  double *logps;
  PyObject *logP;
  
  static char *kwlist[] = {"logP", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &logP))
    return -1;
  
  N = PyArray_DIM(logP, 0);
  K = PyArray_DIM(logP, 0);
  logps = (double *)PyArray_DATA((PyArrayObject *)logP);
  
  printf("N = %d, K = %c\n", N, K);
  
  int i, j;
  for (i = 1; i < N; ++i)
    for (i = 1; i < K; ++i)
      printf("%f", logps[N*i+j]);
  
  Py_RETURN_NONE;
}

static PyMethodDef StateMethods[] = {
    {"view",  view, METH_VARARGS, "View the logP matrix."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef statemodule = {
    PyModuleDef_HEAD_INIT,
    "state",
    spam_doc,
    -1,
    StateMethods
};

PyMODINIT_FUNC
PyInit_state(void)
{
    return PyModule_Create(&statemodule);
}


