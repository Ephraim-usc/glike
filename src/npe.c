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
  double *data;
  PyObject *logP;
  
  static char *kwlist[] = {"logP", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &logP))
    Py_RETURN_NONE;
  
  N = PyArray_DIM(logP, 0);
  K = PyArray_DIM(logP, 1);
  data = (double *)PyArray_DATA((PyArrayObject *)logP);
  
  //computing nums and num
  int num = 1;
  int *nums = calloc(N, sizeof(int));
  for (n = 0; n < N; n++)
  {
    for (k = 0; k < K; k++)
      if(data[K*n+k] > -INFINITY) nums[n]++;
    num *= nums[n];
  }
  
  // computing values
  int *values = (int *)malloc(num * N * sizeof(int)); int *values_;
  int *logps = (double *)malloc(num * N * sizeof(double)); int *logps_;
  int size = num; int chunk;
  for (n = 0; n < N; n++)
  {
    values_ = values + num * n;
    logps_ = logps + num * n;
    p = values_; q = logps_;
    size /= nums[n];
    
    // memset
    int offset = 0; int end;
    double datum;
    for (k = 0; k < K; k++)
    {
      datum = data[K*n+k];
      if(data[K*n+k] == -INFINITY) continue;
      for (end = offset + size; offset < end; offset++)
      {
          values_[offset] = k;
          logps_[offset] = datum;
      }
    }
    
    // memcpy
    chunk = offset;
    for (; offset < num; offset += chunk)
    {
      memcpy(values_ + offset, values_, chunk * sizeof(int));
      memcpy(logps_ + offset, logps_, chunk * sizeof(double));
    }
  }
  
  npy_intp dims[] = {N, num};
  PyObject *values_array = PyArray_SimpleNewFromData(2, dims, NPY_INT, values);
  PyObject *logps_array = PyArray_SimpleNewFromData(2, dims, NPY_INT, logps);
  
  values_array = PyArray_Transpose(values_array, NULL);
  logps_array = PyArray_Transpose(logps_array, NULL);
  
  PyObject *out = PyTuple_Pack(2, values_array, logps_array);
  return out;
}








static PyMethodDef npeMethods[] = {
  {"view", (PyCFunction) view, METH_VARARGS | METH_KEYWORDS, "View the logP matrix."},
  {"product_det", (PyCFunction) product_det, METH_VARARGS | METH_KEYWORDS, "Deterministic product."},
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef npemodule = {
  PyModuleDef_HEAD_INIT,
  "npe",
  "npe Module",
  -1,
  npeMethods
};

PyMODINIT_FUNC
PyInit_npe(void)
{
  import_array();
  return PyModule_Create(&npemodule);
}


