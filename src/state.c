#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "structmember.h"




typedef struct StateObject
{
  PyObject_HEAD
  int len;
  int *values;
  struct StateObject **parents;
  struct StateObject **children;
} StateObject;

static void State_dealloc(StateObject *self)
{
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *State_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  StateObject *self;
  self = (StateObject *) type->tp_alloc(type, 0);
  return (PyObject *) self;
}

static int State_init(StateObject *self, PyObject *args, PyObject *kwds)
{
  self->len = 0;
  self->values = NULL;
  self->parents = NULL;
  self->children = NULL;
  
  PyObject *values;
  
  static char *kwlist[] = {"values", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &values))
    return -1;
  
  self->len = PyArray_DIM(values, 0);
  
  int *p = (int *)PyArray_DATA((PyArrayObject *)values);
  self->values = (int *)malloc(self->len * sizeof(int));
  memcpy(self->values, p, self->len * sizeof(int));
  
  return 0;
}

static PyObject *State_print(StateObject *self, PyObject *args)
{
  int len = self->len;
  
  int i;
  for (i = 0; i < len; i++)
    printf("%d ", self->values[i]);
  printf("\n");
  
  Py_RETURN_NONE;
}

static PyMethodDef State_methods[] = 
{
  {"print", (PyCFunction) State_print, METH_NOARGS, "print state"},
  {NULL},
};

static PyTypeObject StateType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "state.State",
    .tp_basicsize = sizeof(StateObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = State_new,
    .tp_init = (initproc) State_init,
    .tp_dealloc = (destructor) State_dealloc,
    .tp_methods = State_methods,
};





typedef struct
{
  PyObject_HEAD
  double t;
  int dim_in;
  int dim_out;
  double *logP;
  double **logRR;
} TransitionObject;

static void Transition_dealloc(TransitionObject *self)
{
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *Transition_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  TransitionObject *self;
  self = (TransitionObject *) type->tp_alloc(type, 0);
  return (PyObject *) self;
}

static int Transition_init(TransitionObject *self, PyObject *args, PyObject *kwds)
{
  self->t = 0;
  self->dim_in = 0;
  self->dim_out = 0;
  self->logP = NULL;
  self->logRR = NULL;
  
  PyObject *logP;
  PyObject *logRR;
  
  static char *kwlist[] = {"t", "dim_in", "dim_out", "logP", "logRR", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|diiOO", kwlist, &self->t, &self->dim_in, &self->dim_out, &logP, &logRR))
    return -1;
  
  int dim_in = self->dim_in;
  int dim_out = self->dim_out;
  
  if ((PyArray_DIM(logP, 0) != dim_in * dim_out)||
      (PyArray_DIM(logRR, 0) != dim_in * dim_out)||
      (PyArray_DIM(logRR, 1) != dim_in * dim_out))
  {
    printf("error: dimensions do not match!\n");
    return -1;
  }
  
  double *p = (double *)PyArray_DATA((PyArrayObject *)logP);
  self->logP = (double *)malloc(dim_in * dim_out * sizeof(double));
  memcpy(self->logP, p, dim_in * dim_out * sizeof(double));
  
  int i;
  double *rr = (double *)PyArray_DATA((PyArrayObject *)logRR);
  self->logRR = (double **)malloc(dim_in * dim_out * sizeof(double *));
  for (i = 0; i < dim_in * dim_out; i++)
  {
    self->logRR[i] = (double *)malloc(dim_in * dim_out * sizeof(double));
    memcpy(self->logRR[i], rr + i * dim_in * dim_out, dim_in * dim_out * sizeof(double));
  }
  
  return 0;
}

static PyObject *Transition_print(TransitionObject *self, PyObject *args)
{
  int dim_in = self->dim_in;
  int dim_out = self->dim_out;
  printf("%9.4lf\n\n", self->t);
  
  int i;
  int j;
  
  for (i = 0; i < dim_in * dim_out; i++)
    printf("       %d%d", i/dim_in, i%dim_in);
  printf("\n\n");
  
  for (i = 0; i < dim_in * dim_out; i++)
    printf("%9.4lf", self->logP[i]);
  printf("\n\n");
  
  for (i = 0; i < dim_in * dim_out; i++)
  {
    for (j = 0; j < dim_in * dim_out; j++)
      printf("%9.4lf", self->logRR[i][j]);
    printf("\n");
  }
  
  Py_RETURN_NONE;
}

static PyMethodDef Transition_methods[] = 
{
  {"print", (PyCFunction) Transition_print, METH_NOARGS, "print transition"},
  {NULL},
};


static PyTypeObject TransitionType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "state.Transition",
    .tp_basicsize = sizeof(TransitionObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = Transition_new,
    .tp_init = (initproc) Transition_init,
    .tp_dealloc = (destructor) Transition_dealloc,
    .tp_methods = Transition_methods,
};



static struct PyModuleDef stateModule =
{
  PyModuleDef_HEAD_INIT,
  .m_name = "stateModule",
  .m_doc = "state Module",
  .m_size = -1,
};

PyMODINIT_FUNC 
PyInit_state(void)
{
  if (PyType_Ready(&TransitionType) < 0)
    return NULL;
  if (PyType_Ready(&StateType) < 0)
    return NULL;
  
  PyObject *m = PyModule_Create(&stateModule);
  if (m == NULL)
    return NULL;
  
  Py_INCREF(&TransitionType);
  if (PyModule_AddObject(m, "Transition", (PyObject *) &TransitionType) < 0) 
  {
    Py_DECREF(&TransitionType);
    Py_DECREF(m);
    return NULL;
  }
  
  Py_INCREF(&StateType);
  if (PyModule_AddObject(m, "State", (PyObject *) &StateType) < 0) 
  {
    Py_DECREF(&StateType);
    Py_DECREF(m);
    return NULL;
  }
  
  return m;
}

