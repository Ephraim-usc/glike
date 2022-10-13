#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "structmember.h"


typedef struct
{
  PyObject_HEAD
  double t;
} TransitionObject;

static void Transition_dealloc(TransitionObject *self)
{
  Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *Transition_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  TransitionObject *self;
  self = (TransitionObject *) type->tp_alloc(type, 0);
  if (self != NULL) 
  {
    self->t = 0;
  }
  return (PyObject *) self;
}

static int Transition_init(TransitionObject *self, PyObject *args, PyObject *kwds)
{
  self->t = 0;
  
  static char *kwlist[] = {"t", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|d", kwlist, &self->t))
    return -1;
  
  return 0;
}

static PyMemberDef Transition_members[] = 
{
    {"t", T_DOUBLE, offsetof(TransitionObject, t), 0, "transition time"},
    {NULL}  /* Sentinel */
};

static PyObject * Transition_print(TransitionObject *self, PyObject *args)
{
    return printf("%lf\n", self->t);
}

static PyMethodDef Transition_methods[] = 
{
  {"print", (PyCFunction) Transition_print, METH_NOARGS, "print transition"},
  {NULL},
};


static PyTypeObject CustomType = {
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
  if (PyType_Ready(&CustomType) < 0)
    return NULL;
  
  m = PyModule_Create(&statemodule);
  if (m == NULL)
    return NULL;
  
  Py_INCREF(&CustomType);
  if (PyModule_AddObject(m, "Custom", (PyObject *) &CustomType) < 0) 
  {
    Py_DECREF(&CustomType);
    Py_DECREF(m);
    return NULL;
  }
  
  return m;
}

