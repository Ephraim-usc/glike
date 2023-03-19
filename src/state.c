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
    Py_RETURN_NONE;
  
  N = PyArray_DIM(logP, 0);
  K = PyArray_DIM(logP, 1);
  logps = (double *)PyArray_DATA((PyArrayObject *)logP);
  
  printf("N = %d, K = %d\n", N, K);
  
  int i, j;
  for (i = 0; i < N; i++)
    for (j = 0; j < K; j++)
      printf("%f ", logps[K*i+j]);
  printf("\n");
  
  Py_RETURN_NONE;
}

static PyObject *product_det(PyObject *self, PyObject *args, PyObject *kwds)
{
  int N,K;
  int n, k;
  double *logps;
  PyObject *logP;
  
  static char *kwlist[] = {"logP", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &logP))
    Py_RETURN_NONE;
  
  N = PyArray_DIM(logP, 0);
  K = PyArray_DIM(logP, 1);
  logps = (double *)PyArray_DATA((PyArrayObject *)logP);
  
  printf("N = %d, K = %d\n", N, K);
  
  int num = 1, num_ = 0;
  int *nums = calloc(N, sizeof(int));
  for (n = 0; n < N; n++)
  {
    for (k = 0; k < K; k++)
      if(logps[K*n+k] > -INFINITY) nums[n]++;
    num *= nums[n];
  }
  printf("num = %d\n", num);
  
  int *p, *q;
  int *values = (int *)malloc(num * N * sizeof(int));
  int size = num;
  for (n = 0; n < N; n++)
  {
    values_ = values + num * n;
    p = values_;
    
    for (k = 0; k < K; k++)
      if(logps[K*n+k] > -INFINITY)
        for (q = p + size/nums[n]; p < q; p++)
          *p = k;
    
    for (q = values_ + num; p < q; p += size)
      memcpy(p, values_, size * sizeof(int))
    
    size /= nums[n];
  }
  
  
  Py_RETURN_NONE;
}








static PyMethodDef StateMethods[] = {
    {"view", (PyCFunction) view, METH_VARARGS | METH_KEYWORDS, "View the logP matrix."},
    {"product_det", (PyCFunction) product_det, METH_VARARGS | METH_KEYWORDS, "Deterministic product."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef statemodule = {
    PyModuleDef_HEAD_INIT,
    "state",
    "state Module",
    -1,
    StateMethods
};

PyMODINIT_FUNC
PyInit_state(void)
{
    return PyModule_Create(&statemodule);
}


