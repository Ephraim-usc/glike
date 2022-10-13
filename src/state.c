#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "structmember.h"


typedef struct transition 
{
  double t;
} transition;

transition* new_transition(double t)
{
  transition* trn;
  trn = (transition *)malloc(sizeof(transition));
  trn->t = t;
  return trn;
}

void print_transition(transition *trn)
{
  printf("%lf", trn->t);
}



static PyObject* py_new_transition(PyObject* self, PyObject* args)
{ 
  double t;
  PyArg_ParseTuple(args, "d", &t);
  
  transition* trn = new_transition((int)t);
  PyObject* py_trn = PyCapsule_New((void *)trn, "state._transition_C_API", NULL);
  
  return py_trn;
}

static PyObject* py_print_transition(PyObject* self, PyObject* args)
{
  PyObject* py_trn;
  PyArg_UnpackTuple(self, NULL, 1, 1, &py_trn);
  transition* trn = (transition *)PyCapsule_GetPointer(py_trn, "state._transition_C_API");
  print_transition(trn);
  
  Py_RETURN_NONE;
}




static PyMethodDef myMethods[] = 
{
  {"new_transition", py_new_transition, METH_VARARGS, "new transition"},
  {"print", py_print_transition, METH_VARARGS, "print transition"},
  {NULL, NULL, 0, NULL},
};

static struct PyModuleDef stateModule =
{
  PyModuleDef_HEAD_INIT,
  "stateModule",
  "state Module",
  -1,
  myMethods
};

PyMODINIT_FUNC PyInit_state(void)
{
  return PyModule_Create(&stateModule);
}

