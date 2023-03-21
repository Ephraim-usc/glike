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
  int N, K;
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
  double *logps = (double *)malloc(num * N * sizeof(double)); double *logps_;
  int size = num; int chunk;
  for (n = 0; n < N; n++)
  {
    values_ = values + num * n;
    logps_ = logps + num * n;
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
  PyObject *logps_array = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, logps);
  
  values_array = PyArray_Transpose((PyArrayObject *)values_array, NULL);
  logps_array = PyArray_Transpose((PyArrayObject *)logps_array, NULL);
  
  PyObject *out = PyTuple_Pack(2, values_array, logps_array);
  return out;
}

static PyObject *product_rand(PyObject *self, PyObject *args, PyObject *kwds)
{
  int N, K, m;
  int n, k;
  double *data, *data_;
  PyObject *P;
  
  static char *kwlist[] = {"P", "", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Oi", kwlist, &P, &m))
    Py_RETURN_NONE;
  
  N = PyArray_DIM(P, 0);
  K = PyArray_DIM(P, 1);
  data = (double *)PyArray_DATA((PyArrayObject *)P);
  
  double *cdf = (double *)malloc(N * K * sizeof(double)); double *cdf_;
  int *values = (int *)malloc(N * K * sizeof(int)); int *values_;
  
  int i, j;
  for (n = 0; n < N; n++)
  {
    data_ = data + n * K;
    cdf_ = cdf + n * K;
    values_ = values + n * K;
    cdf_[0] = 0.0;
    i = 0;
    for (k = 0; k < K; k++)
    {
      if (data_[k] <= 0.0) continue;
      cdf_[i] = data_[k];
      values_[i] = k;
      i ++;
    }
    for (j = 1; j < i; j++)
    {
      cdf_[j] += cdf_[j-1];
    }
  }
  
  for (i = 0; i < N*K; i++)
    printf("%f ", cdf[i]);
  printf("\n");
  for (i = 0; i < N*K; i++)
    printf("%d ", values[i]);
  printf("\n");
  
  Py_RETURN_NONE;
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


