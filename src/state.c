#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "structmember.h"


typedef struct
{
  PyObject_HEAD
  double t;
  int dim_in;
  int dim_out;
  double *logP;
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
  
  PyObject *logP;
  
  static char *kwlist[] = {"t", "logP", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dO", kwlist, &self->t, &logP))
    return -1;
  
  double *p = (double *)PyArray_DATA((PyArrayObject *)logP);
  self->dim_in = PyArray_DIM(logP, 0);
  self->dim_out = PyArray_DIM(logP, 1);
  
  self->logP = (double *)malloc(self->dim_in * self->dim_out * sizeof(double));
  memcpy(self->logP, p, self->dim_in * self->dim_out * sizeof(double));
  
  return 0;
}

static PyMemberDef Transition_members[] = 
{
  {"t", T_DOUBLE, offsetof(TransitionObject, t), 0, "transition time"},
  {NULL}  /* Sentinel */
};

static PyObject *Transition_print(TransitionObject *self, PyObject *args)
{
  printf("%lf\n", self->t);
  
  int i;
  int j;
  for (i = 0; i < self->dim_in; i++)
  {
    for (j = 0; j < self->dim_out; j++)
      printf("%lf", self->logP[i * self->dim_in + j]);
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
    .tp_members = Transition_members,
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
  PyObject *m;
  if (PyType_Ready(&TransitionType) < 0)
    return NULL;
  
  m = PyModule_Create(&stateModule);
  if (m == NULL)
    return NULL;
  
  Py_INCREF(&TransitionType);
  if (PyModule_AddObject(m, "Transition", (PyObject *) &TransitionType) < 0) 
  {
    Py_DECREF(&TransitionType);
    Py_DECREF(m);
    return NULL;
  }
  
  return m;
}

