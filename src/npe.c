#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

void free_wrap(PyObject *capsule) {
    void *memory = PyCapsule_GetPointer(capsule, NULL);
    free(memory);
}

static PyObject *free_(PyObject *self, PyObject *args, PyObject *kwds)
{
  PyObject *x;
  void *data;
  
  static char *kwlist[] = {"x", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &x))
    Py_RETURN_NONE;
  
  data = (void *)PyArray_DATA((PyArrayObject *)x);
  free(data);
  
  Py_RETURN_NONE;
}


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
  
  // make buffers
  npy_intp dims[] = {num, N};
  
  PyObject *values_array = PyArray_Zeros(2, dims, PyArray_DescrFromType(NPY_INT), 0);
  PyObject *logps_array = PyArray_Zeros(2, dims, PyArray_DescrFromType(NPY_DOUBLE), 0);
  
  PyArray_STRIDES(values_array)[0] = sizeof(int); PyArray_STRIDES(values_array)[1] = num * sizeof(int); 
  PyArray_STRIDES(logps_array)[0] = sizeof(double); PyArray_STRIDES(logps_array)[1] = num * sizeof(double); 
  
  int *values = (int *)PyArray_DATA((PyArrayObject *)values_array); int *values_;
  double *logps = (double *)PyArray_DATA((PyArrayObject *)logps_array); double *logps_;
  
  // computing values
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
  
  PyObject *out = PyTuple_Pack(2, values_array, logps_array);
  
  // this is required since PyTuple_Pack increments ref count of each element.
  Py_DECREF(values_array);
  Py_DECREF(logps_array);
  return out;
}

// it's required that P is row-first allocated in memory
static PyObject *product_sto(PyObject *self, PyObject *args, PyObject *kwds)
{
  int N, K, M;
  int n, k, m;
  double *data, *data_;
  PyObject *P, *values_array, *ps_array;
  
  static char *kwlist[] = {"P", "num", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OiOO", kwlist, &P, &M, &values_array, &ps_array))
    Py_RETURN_NONE;
  
  N = PyArray_DIM(P, 0);
  K = PyArray_DIM(P, 1);
  data = (double *)PyArray_DATA((PyArrayObject *)P);
  
  // preparations
  double *pdf = (double *)malloc(N * K * sizeof(double)); double *pdf_;
  double *cdf = (double *)malloc(N * K * sizeof(double)); double *cdf_;
  int *idx = (int *)malloc(N * K * sizeof(int)); int *idx_;
  
  int i, j;
  for (n = 0; n < N; n++)
  {
    data_ = data + n * K;
    pdf_ = pdf + n * K;
    cdf_ = cdf + n * K;
    idx_ = idx + n * K;
    
    i = 0;
    for (k = 0; k < K; k++)
    {
      if (data_[k] < 1e-8) continue;
      pdf_[i] = data_[k];
      idx_[i] = k;
      i ++;
    }
    
    cdf_[0] = pdf_[0];
    for (j = 1; j < i; j++)
    {
      cdf_[j] = cdf_[j-1] + pdf_[j];
    }
    
    if (fabs(cdf_[j-1] - 1.0) > 1e-8)
    {
      printf("Error: probabilities don't sum up to 1!\n");
      Py_RETURN_NONE;
    }
  }
  
  /*
  for (i = 0; i < N*K; i++)
    printf("%f ", pdf[i]);
  printf("\n");
  for (i = 0; i < N*K; i++)
    printf("%f ", cdf[i]);
  printf("\n");
  for (i = 0; i < N*K; i++)
    printf("%d ", idx[i]);
  printf("\n"); 
  */
  
  // making buffers
  PyArray_DIMS(values_array)[0] = M; PyArray_DIMS(values_array)[1] = N; 
  PyArray_DIMS(ps_array)[0] = M; PyArray_DIMS(ps_array)[1] = N; 
  
  PyArray_STRIDES(values_array)[0] = sizeof(int); PyArray_STRIDES(values_array)[1] = M * sizeof(int); 
  PyArray_STRIDES(ps_array)[0] = sizeof(double); PyArray_STRIDES(ps_array)[1] = M * sizeof(double); 
  
  int *values = (int *)PyArray_DATA((PyArrayObject *)values_array); int *values_;
  double *ps = (double *)PyArray_DATA((PyArrayObject *)ps_array); double *ps_;  
  
  // computing values
  double tmp;
  for (n = 0; n < N; n++)
  {
    data_ = data + n * K;
    pdf_ = pdf + n * K;
    cdf_ = cdf + n * K;
    idx_ = idx + n * K;
    values_ = values + n * M;
    ps_ = ps + n * M;
    
    for (m = 0; m < M; m++)
    {
      tmp = drand48();
      i = 0; while(cdf_[i] < tmp) i++;
      values_[m] = idx_[i];
      ps_[m] = pdf_[i];
    }
  }
  
  free(pdf);
  free(cdf);
  free(idx);
  
  PyObject *out = PyTuple_Pack(2, values_array, ps_array);
  
  // this is required since PyTuple_Pack increments ref count of each element.
  Py_DECREF(values_array); 
  Py_DECREF(ps_array);
  return out;
}




static PyMethodDef npeMethods[] = {
  {"view", (PyCFunction) view, METH_VARARGS | METH_KEYWORDS, "View the logP matrix."},
  {"product_det", (PyCFunction) product_det, METH_VARARGS | METH_KEYWORDS, "Deterministic product."},
  {"product_sto", (PyCFunction) product_sto, METH_VARARGS | METH_KEYWORDS, "Stochastic product."},
  {"free", (PyCFunction) free_, METH_VARARGS | METH_KEYWORDS, "Manually free an array."},
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


