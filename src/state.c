#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "list.h"

typedef double DTYPE;
typedef unsigned long ITYPE;










typedef struct state
{
  ITYPE *values; // population labels of the lineages
  DTYPE logp;
  
  list *logps_anc;
  plist *links_anc;
  
  list *logps_des;
  plist *links_des;
} state;

state* new_state(ITYPE *values)
{
  state* stt;
  stt = (state *)malloc(sizeof(state));
  stt->values = values;
  stt->logp = 1.0; // any logp greater than 0 means nan.
  
  stt->num_anc = 0;
  stt->num_des = 0;
  return stt;
}

void print_state(state* stt)
{
  ITYPE n = stt->n;
  DTYPE t = stt->t;
  ITYPE *values = stt->values;
  
  printf("state of %lu lineages at time %lf\n", n, t);
  
  int i;
  for(i = 0; i < n; i++)
    printf("%lu", *(values+i));
  printf("\n");
  
  printf("ancestral links: ")
  for(i = 0; i < stt->num_anc; i++)
    printf("%lf", *(stt->logps_anc+i));
  printf("\n");
  
  printf("descendent links: ")
  for(i = 0; i < stt->num_des; i++)
    printf("%lf", *(stt->logps_des+i));
  printf("\n");
}

typedef struct bundle
{
  ITYPE n; // number of lineages of the state
  DTYPE t; // time of the state
  ITYPE *lins; // lineages of the state
  state **states;
} bundle;

bundle* new_bundle(ITYPE n, DTYPE t, ITYPE *lins)
{
  bundle* bdl;
  bdl = (bundle *)malloc(sizeof(bundle));
  bdl->n = n;
  bdl->t = t;
  bdl->lins = lins;
  
  return bdl;
}



static PyMethodDef myMethods[] = 
{
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

/*
int main()
{
  ITYPE lins[5] = {4,9,2,10,31};
  ITYPE values[5] = {0,1,1,2,0};
  state *stt = new_state(5, 37.8, lins, values);
  print_state(stt);
  
  return 0;
}
*/
