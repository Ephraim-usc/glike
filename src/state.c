#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "structmember.h"













typedef struct State
{
  PyObject_HEAD
  int len;
  int *values;
  struct State **parents;
  struct State **children;
} State;

State *State_new()
{
  State *state = (State *)malloc(sizeof(State));
  state->len = 0;
  state->values = NULL;
  state->parents = NULL;
  state->children = NULL;
  return state;
}







typedef struct BundleObject
{
  PyObject_HEAD
  double t;
  int len;
  int num_states;
  int *lineages;
  State **states;
} BundleObject;

static void Bundle_dealloc(BundleObject *self)
{
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *Bundle_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  BundleObject *self;
  self = (BundleObject *) type->tp_alloc(type, 0);
  return (PyObject *) self;
}

static int Bundle_init(BundleObject *self, PyObject *args, PyObject *kwds)
{
  self->t = 0;
  self->len = 0;
  self->lineages = NULL;
  self->states = NULL;
  
  PyObject *lineages;
  PyObject *values;
  
  static char *kwlist[] = {"t", "lineages", "values", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dOO", kwlist, &self->t, &lineages, &values))
    return -1;
  
  int len = PyArray_DIM(lineages, 0);
  int num_states = PyArray_DIM(values, 0);
  self->len = len;
  self->num_states = num_states;
  
  if (len != PyArray_DIM(values, 1))
  {
    printf("error: dimensions do not match!\n");
    return -1;
  }
  
  int *l = (int *)PyArray_DATA((PyArrayObject *)lineages);
  self->lineages = (int *)malloc(len * sizeof(int));
  memcpy(self->lineages, l, len * sizeof(int));
  
  int i; State *state;
  int *v = (int *)PyArray_DATA((PyArrayObject *)values);
  self->states = (State **)malloc(num_states * sizeof(State *));
  for (i = 0; i < num_states; i++)
  {
    state = State_new();
    state->len = len;
    state->values = (int *)malloc(len * sizeof(int));
    memcpy(state->values, v + i * len, len * sizeof(int));
    self->states[i] = state;
  }
  
  return 0;
}


static PyObject *Bundle_print(BundleObject *self, PyObject *args)
{
  double t = self->t;
  int len = self->len;
  int num_states = self->num_states;
  
  printf("%lf\n\n", self->t);
  
  int i, j;
  
  for (i = 0; i < len; i++)
    printf("%d ", self->lineages[i]);
  printf("\n\n");
  
  for (i = 0; i < num_states; i++)
  {
    for (j = 0; j < len; j++)
      printf("%d ", self->states[i]->values[j]);
    printf("\n");
  }
  
  Py_RETURN_NONE;
}

static PyMethodDef Bundle_methods[] = 
{
  {"print", (PyCFunction) Bundle_print, METH_NOARGS, "print state bundle"},
  {NULL},
};


static PyTypeObject BundleType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "state.Bundle",
    .tp_basicsize = sizeof(BundleObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_new = Bundle_new,
    .tp_init = (initproc) Bundle_init,
    .tp_dealloc = (destructor) Bundle_dealloc,
    .tp_methods = Bundle_methods,
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
  
  int i, j;
  
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
  if (PyType_Ready(&BundleType) < 0)
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
  
  Py_INCREF(&BundleType);
  if (PyModule_AddObject(m, "Bundle", (PyObject *) &BundleType) < 0) 
  {
    Py_DECREF(&BundleType);
    Py_DECREF(m);
    return NULL;
  }
  
  return m;
}

