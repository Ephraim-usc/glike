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
  
  //computing nums and num
  int num = 1, num_ = 0;
  int *nums = calloc(N, sizeof(int));
  for (n = 0; n < N; n++)
  {
    for (k = 0; k < K; k++)
      if(logps[K*n+k] > -INFINITY) nums[n]++;
    num *= nums[n];
  }
  
  // computing values
  int *p, *q;
  int *values = (int *)malloc(num * N * sizeof(int)); int *values_;
  int size = num; int chunk;
  for (n = 0; n < N; n++)
  {
    values_ = values + num * n;
    p = values_;
    size /= nums[n];
    
    // memset
    for (k = 0; k < K; k++)
      if(logps[K*n+k] > -INFINITY)
        for (q = p + size; p < q; p++)
          *p = k;
    
    // memcpy
    chunk = size * nums[n];
    for (q = values_ + num; p < q; p += chunk)
      memcpy(p, values_, chunk * sizeof(int));
  }
  
  
  int i, j;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < num; j++)
      printf("%d ", values[i*num + j]);
    printf("\n");
  }
  
  
  int dims[2]; dims[0] = 2; dims[1] = 3; 
  int data[6] = {1,2,3,4,5,6};
  Py_Initialize();
  import_array();
  PyObject *out = PyArray_SimpleNew(2, dims, NPY_FLOAT);
  
  return out;
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


