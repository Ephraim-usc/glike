#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "structmember.h"

#include "list.h"

static double dot(double *x, int *y, int len) // note that y is (int *), this is just for our algorithm.
{
  double buffer = 0;
  int i;
  for (i = 0; i < len; i++)
  {
    if (isnan(x[i]))
      continue;
    buffer += x[i] * y[i];
  }
  return buffer;
}









typedef struct State
{
  PyObject_HEAD
  int len;
  double logp;
  int *values;
  int num_parents;
  int num_children;
  double *logps_parents;
  double *logps_children;
  struct State **parents;
  struct State **children;
} State;

State *State_new()
{
  State *state = (State *)malloc(sizeof(State));
  state->len = 0;
  state->logp = -INFINITY;
  state->values = NULL;
  state->num_parents = 0;
  state->num_children = 0;
  state->logps_parents = NULL;
  state->logps_children = NULL;
  state->parents = NULL;
  state->children = NULL;
  return state;
}

void State_free(State *state)
{
  if (state->values) free(state->values);
  if (state->logps_parents) free(state->logps_parents);
  if (state->logps_children) free(state->logps_children);
  if (state->parents) free(state->parents);
  if (state->children) free(state->children);
  free(state);
}

void State_print(State *state)
{
  printf("[%d parents]  ", state->num_parents);
  
  printf("logp=%.4lf ", state->logp);
  printf("< ");
  int i;
  for (i = 0; i < state->len; i++)
    printf("%d ", state->values[i]);
  printf(">");
  
  printf("  [%d children]\n", state->num_children);
}




typedef struct BundleObject
{
  PyObject_HEAD
  double t;
  int len;
  int num_states;
  int *lineages;
  State **states;
  struct BundleObject *parent;
  struct BundleObject *child;
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
    state->logp = -INFINITY;
    state->values = (int *)malloc(len * sizeof(int));
    memcpy(state->values, v + i * len, len * sizeof(int));
    self->states[i] = state;
  }
  
  return 0;
}

static PyObject *Bundle_free(BundleObject *self)
{
  BundleObject *bundle = self;
  while (bundle != NULL)
  {
    int i;
    for (i = 0; i < bundle->num_states; i++)
    {
      State_free(bundle->states[i]);
    }
    if (bundle->lineages) free(bundle->lineages);
    if (bundle->states) free(bundle->states);
    
    bundle = bundle->child;
  }
  
  Py_RETURN_NONE;
}

static PyObject *Bundle_print(BundleObject *self, PyObject *args)
{
  int len = self->len;
  int num_states = self->num_states;
  
  printf("%lf\n\n", self->t);
  
  int i;
  for (i = 0; i < len; i++)
    printf("%d ", self->lineages[i]);
  printf("\n\n");
  
  for (i = 0; i < num_states; i++)
  {
    State_print(self->states[i]);
  }
  
  Py_RETURN_NONE;
}

static PyObject *Bundle_propagate(BundleObject *self, PyObject *args)
{
  int i; // index of state
  int j; // index of link
  for (i = 0; i < self->num_states; i++)
    self->states[i]->logp = 0;
  
  BundleObject *bundle = self->parent;
  while (bundle != NULL)
  {
    for (i = 0; i < bundle->num_states; i++)
    {
      State *state = bundle->states[i];
      if (state->num_children == 0)
        continue;
      
      double *buffer = (double *)malloc(state->num_children * sizeof(double));
      memcpy(buffer, state->logps_children, state->num_children * sizeof(double));
      
      double tmp = 0;
      double current_max = -INFINITY;
      for (j = 0; j < state->num_children; j++)
      {
        *(buffer + j) += state->children[j]->logp;
        if (*(buffer + j) > current_max) current_max = *(buffer + j);
      }
      
      for (j = 0; j < state->num_children; j++)
      {
        tmp += exp(*(buffer + j) - current_max);
      }
      
      state->logp = current_max + log(tmp);
      free(buffer);
    }
    
    bundle = bundle->parent;
  }
  
  Py_RETURN_NONE;
}

static PyObject *Bundle_logp(BundleObject *self, PyObject *args)
{
  int s;
  double tmp = 0;
  double current_max = -INFINITY;
  double logp_;
  for (s = 0; s < self->num_states; s++)
  {
    logp_ = self->states[s]->logp;
    if (logp_ > current_max) current_max = logp_;
  }
  for (s = 0; s < self->num_states; s++)
  {
    logp_ = self->states[s]->logp;
    tmp += exp(logp_ - current_max);
  }
  
  return Py_BuildValue("d", current_max + log(tmp));
}


static PyObject *Bundle_diverge(BundleObject *self, PyObject *args, PyObject *kwds);

static PyObject *Bundle_evolve(BundleObject *self, PyObject *args, PyObject *kwds); // inverse of transition

static PyMethodDef Bundle_methods[] = 
{
  {"free", (PyCFunction) Bundle_free, METH_NOARGS, "release memory of this and all descendent bundles"},
  {"print", (PyCFunction) Bundle_print, METH_NOARGS, "print state bundle"},
  {"propagate", (PyCFunction) Bundle_propagate, METH_NOARGS, "logp propagation from this bundle"},
  {"logp", (PyCFunction) Bundle_logp, METH_NOARGS, "compute the total logp of this bundle"},
  {"diverge", (PyCFunction) Bundle_diverge, METH_VARARGS | METH_KEYWORDS, "bundle diverge"},
  {"evolve", (PyCFunction) Bundle_evolve, METH_VARARGS | METH_KEYWORDS, "bundle evolve (inverse of transition)"},
  {NULL},
};


static PyObject *Bundle_getparent(BundleObject *self, void *closure)
{
  if (self->parent)
  {
    Py_INCREF(self->parent);
    return (PyObject *)self->parent;
  }
  else
    Py_RETURN_NONE;
}

static PyObject *Bundle_getchild(BundleObject *self, void *closure)
{
  if (self->child)
  {
    Py_INCREF(self->child);
    return (PyObject *)self->child;
  }
  else
    Py_RETURN_NONE;
}

static PyGetSetDef Bundle_getsetters[] = {
    {"parent", (getter) Bundle_getparent, NULL, "parent bundle", NULL},
    {"child", (getter) Bundle_getchild, NULL, "child bundle", NULL},
    {NULL}  /* Sentinel */
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
    .tp_getset = Bundle_getsetters,
};


static PyObject *Bundle_diverge(BundleObject *self, PyObject *args, PyObject *kwds)
{
  PyObject *parent;
  PyObject *children;
  PyObject *logn;
  
  static char *kwlist[] = {"parent", "children", "logn", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OOO", kwlist, &parent, &children, &logn))
    Py_RETURN_NONE;
  
  if (PyArray_DIM(parent, 0) != 1)
  {
    printf("error: parent must be of length 1!\n");
    Py_RETURN_NONE;
  }
  
  int num_children = PyArray_DIM(children, 0);
  int *p = (int *)PyArray_DATA((PyArrayObject *)parent);
  int *c = (int *)PyArray_DATA((PyArrayObject *)children);
  double *l = (double *)PyArray_DATA((PyArrayObject *)logn);
  
  // find index of parent
  int i;
  int index = -1;
  for (i = 0; i<self->len; i++)
  {
    if (self->lineages[i] == *p)
    {
      index = i;
      break;
    }
  }
  
  if (index == -1)
  {
    printf("error: parent lineage not found!\n");
    Py_RETURN_NONE;
  }
  
  // create new bundle with lineages (children lineages are added at the end)
  BundleObject *bundle = (BundleObject *) BundleType.tp_alloc(&BundleType, 0);
  bundle->t = self->t;
  bundle->len = self->len + num_children - 1;
  bundle->num_states = self->num_states;
  
  bundle->lineages = (int *)malloc((bundle->len) * sizeof(int));
  memcpy(bundle->lineages, self->lineages, index * sizeof(int));
  memcpy(bundle->lineages + index + num_children, self->lineages + index + 1, (self->len - index - 1) * sizeof(int));
  memcpy(bundle->lineages + index, c, num_children * sizeof(int));
  
  // fill out states in new bundle
  int s;
  bundle->states = (State **)malloc(self->num_states * sizeof(State *));
  for (s = 0; s < self->num_states; s++)
  {
    State *state = State_new();
    state->len = self->len + num_children - 1;
    
    int value = self->states[s]->values[index];
    state->values = (int *)malloc(state->len * sizeof(int));
    memcpy(state->values, self->states[s]->values, index * sizeof(int));
    memcpy(state->values + index + num_children, self->states[s]->values + index + 1, (self->len - index - 1) * sizeof(int));
    for (i = index; i < index + num_children; i++)
      state->values[i] = value;
    
    state->num_parents = 1;
    state->parents = (State **)malloc(sizeof(State *));
    state->parents[0] = self->states[s];
    bundle->states[s] = state;
    
    self->states[s]->num_children = 1;
    self->states[s]->children = (State **)malloc(sizeof(State *));
    self->states[s]->children[0] = state;
    self->states[s]->logps_children = (double *)malloc(sizeof(double *));
    self->states[s]->logps_children[0] = l[value];
  }
  
  self->child = bundle;
  bundle->parent = self;
  
  Py_INCREF(bundle);
  return (PyObject *) bundle;
}









typedef struct
{
  PyObject_HEAD
  double t;
  int dim_in;
  int dim_out;
  double *logP;
  double **logRR;
  int *num_outs; // fast track of number of out values for each in value
  int **outs;  // fast track of out values for each in value
  int *num_ins; // fast track of number of in values for each out value
  int **ins;  // fast track of in values for each out value
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
  
  self->num_outs = NULL;
  self->outs = NULL;
  self->num_ins = NULL;
  self->ins = NULL;
  
  PyObject *logP;
  PyObject *logRR;
  
  static char *kwlist[] = {"t", "logP", "logRR", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dOO", kwlist, &self->t, &logP, &logRR))
    return -1;
  
  int dim_in = self->dim_in = PyArray_DIM(logP, 0);
  int dim_out = self->dim_out = PyArray_DIM(logP, 1);
  
  if ((PyArray_DIM(logRR, 0) != dim_in * dim_in)||
      (PyArray_DIM(logRR, 1) != dim_out * dim_out))
  {
    printf("error: dimensions do not match!\n");
    return -1;
  }
  
  double *p = (double *)PyArray_DATA((PyArrayObject *)logP);
  self->logP = (double *)malloc(dim_in * dim_out * sizeof(double));
  memcpy(self->logP, p, dim_in * dim_out * sizeof(double));
  
  int in1, in2, out1, out2;
  p = (double *)PyArray_DATA((PyArrayObject *)logRR);
  self->logRR = (double **)malloc(dim_in * dim_out * sizeof(double *));
  for (in1 = 0; in1 < dim_in; in1++)
  {
    for (out1 = 0; out1 < dim_out; out1++)
    {
      self->logRR[in1 * dim_out + out1] = (double *)malloc(dim_in * dim_out * sizeof(double));
      for (in2 = 0; in2 < dim_in; in2++)
        for (out2 = 0; out2 < dim_out; out2++)
          self->logRR[in1 * dim_out + out1][in2 * dim_out + out2] = p[(in1 * dim_in + in2) * (dim_out * dim_out) + out1 * dim_out + out2];
    }
  }
  
  int in, out, k;
  int *num_outs = (int *)calloc(dim_in, sizeof(int));
  int *num_ins = (int *)calloc(dim_out, sizeof(int));
  int **outs = (int **)malloc(dim_in * sizeof(int *));
  int **ins = (int **)malloc(dim_out * sizeof(int *));
  
  for (in = 0; in < dim_in; in++)
  {
    outs[in] = (int *)malloc(dim_out * sizeof(int));
    k = 0;
    for (out = 0; out < dim_out; out++)
    {
      if (self->logP[in * dim_out + out] == -INFINITY)
        continue;
      num_outs[in] += 1;
      outs[in][k++] = out;
    }
  }
  
  for (out = 0; out < dim_out; out++)
  {
    ins[out] = (int *)malloc(dim_in * sizeof(int));
    k = 0;
    for (in = 0; in < dim_in; in++)
    {
      if (self->logP[in * dim_out + out] == -INFINITY)
        continue;
      num_ins[out] += 1;
      ins[out][k++] = in;
    }
  }
  
  self->num_outs = num_outs;
  self->outs = outs;
  self->num_ins = num_ins;
  self->ins = ins;
  
  return 0;
}

static PyObject *Transition_free(TransitionObject *self)
{
  int i;
  free(self->logP);
  free(self->num_outs);
  free(self->num_ins);
  
  int dim = self->dim_in * self->dim_out;
  for (i = 0; i < dim; i++)
    free(self->logRR[i]);
  free(self->logRR);
  
  for (i = 0; i < self->dim_in; i++)
    free(self->outs[i]);
  free(self->outs);
  
  for (i = 0; i < self->dim_out; i++)
    free(self->ins[i]);
  free(self->ins);
  
  Py_RETURN_NONE;
}

static PyObject *Transition_print(TransitionObject *self, PyObject *args)
{
  int dim_in = self->dim_in;
  int dim_out = self->dim_out;
  printf("%9.4lf\n\n", self->t);
  
  int i, j;
  
  for (i = 0; i < dim_in * dim_out; i++)
    printf("       %d%d", i/dim_out, i%dim_out);
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
  
  /*
  int in, k;
  for (in = 0; in < dim_in; in ++)
  {
    printf("%d ", self->num_outs[in]);
  }
  printf("\n");
  
  for (in = 0; in < dim_in; in ++)
  {
    for (k = 0; k < self->num_outs[in]; k++)
    {
      printf("%d ", self->outs[in][k]);
    }
    printf("\n");
  }
  
  int out, k;
  for (out = 0; out < dim_out; out++)
  {
    printf("%d ", self->num_ins[out]);
  }
  printf("\n");
  
  for (out = 0; out < dim_out; out++)
  {
    for (k = 0; k < self->num_ins[out]; k++)
    {
      printf("%d ", self->ins[out][k]);
    }
    printf("\n");
  }
  */
  
  Py_RETURN_NONE;
}


static PyMethodDef Transition_methods[] = 
{
  {"free", (PyCFunction) Transition_free, METH_NOARGS, "free memory of this transition"},
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


static PyObject *Bundle_evolve(BundleObject *self, PyObject *args, PyObject *kwds)
{
  TransitionObject *transition;
  
  static char *kwlist[] = {"transition", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!", kwlist, &TransitionType, &transition))
    Py_RETURN_NONE;
  
  int len = self->len;
  
  int dim_in = transition->dim_in;
  int dim_out = transition->dim_out;
  int dim = dim_in * dim_out;
  double *logP =  transition->logP;
  double **logRR =  transition->logRR;
  int *num_ins = transition->num_ins;
  int **ins = transition->ins;
  
  Htable *htable = Htable_new(1000);
  
  State *state;
  int s;
  for (s = 0; s < self->num_states; s++)
  {
    state = self->states[s];
    int *values = state->values;
    
    int in, out;
    int i; // index of lineage
    int k; // index of out value
    int z; // position to write
    
    // total number of children states
    int num_children = 1;
    for (i = 0; i < len; i++)
      num_children = num_children * num_ins[values[i]];
    
    // computing three arrays using recursive memory operations
    int *valueses = (int *)calloc(len * num_children, sizeof(int));
    int *Cs = (int *)calloc(dim * num_children, sizeof(int));
    double *logps = (double *)calloc(num_children, sizeof(double));
    
    int current_size_valueses = len;
    int current_size_Cs = dim;
    int current_size = 1;
    for (i = len - 1; i >= 0; i--)
    {
      out = values[i];
      
      for (k = 1; k < num_ins[out]; k++)
      {
        memcpy(valueses + current_size_valueses * k, valueses, current_size_valueses * sizeof(int));
        memcpy(Cs + current_size_Cs * k, Cs, current_size_Cs * sizeof(int));
        memcpy(logps + current_size * k, logps, current_size * sizeof(double));
      }
      
      for (k = 0; k < num_ins[out]; k++)
      {
        in = ins[out][k];
        for (z = current_size * k; z < current_size * (k + 1); z++)
        {
          logps[z] += logP[in * dim_out + out]; // migration probability
          logps[z] += dot(logRR[in * dim_out + out], Cs + z * dim, dim); // non-coalescence probability
        }
        for (z = current_size_valueses * k + i; z < current_size_valueses * (k + 1); z += len)
          valueses[z] = in;
        for (z = current_size_Cs * k + in * dim_out + out; z < current_size_Cs * (k + 1); z += dim)
          Cs[z] += 1;
      }
      current_size_valueses = current_size_valueses * num_ins[out];
      current_size_Cs = current_size_Cs * num_ins[out];
      current_size = current_size * num_ins[out];
    }
    
    // export results from the three memory arrays
    state->num_children = num_children;
    state->logps_children = logps;
    state->children = (State **)malloc(num_children * sizeof(State *));
    
    Hnode *node;
    State *s;
    for (z = 0; z < num_children; z++)
    {
      node = Htable_insert(htable, valueses + len * z, len);
      s = node->pointer;
      
      if (s == NULL)
      {
        s = State_new();
        s->len = len;
        s->values = (int *)malloc(len * sizeof(int));
        memcpy(s->values, valueses + len * z, len * sizeof(int));
        node->pointer = s;
      }
      
      s->num_parents += 1;
      state->children[z] = s;
    }
    
    /*
    for (i = 0; i < len * num_children; i++)
    {
      if (i % len == 0)
        printf(" ");
      printf("%d", valueses[i]);
    }
    printf("\n\n");
    
    for (i = 0; i < dim * num_children; i++)
    {
      if (i % dim == 0)
        printf(" ");
      printf("%d", Cs[i]);
    }
    printf("\n\n");
    
    for (i = 0; i < num_children; i++)
    {
      printf("%lf ", logps[i]);
    }
    printf("\n\n\n");
    */
    
    free(valueses);
    free(Cs);
  }
  
  // creating child bundle
  BundleObject *bundle = (BundleObject *) BundleType.tp_alloc(&BundleType, 0);
  bundle->t = self->t - transition->t;
  bundle->len = self->len;
  bundle->lineages = (int *)malloc(len * sizeof(int));
  memcpy(bundle->lineages, self->lineages, len * sizeof(int));
  
  bundle->states = (State **)Htable_export(htable, &bundle->num_states);
  Htable_free(htable); htable = NULL;
  
  self->child = bundle;
  bundle->parent = self;
  
  Py_INCREF(bundle);
  return (PyObject *) bundle;
}



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

