#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "numpy/arrayobject.h"

typedef double DTYPE;
typedef unsigned long ITYPE;

typedef struct state
{
  ITYPE n; // number of lineages of the state
  DTYPE t; // time of the state
  ITYPE *lins; // lineages of the state
  DTYPE *values; // population labels of the lineages
  
  void *links_anc; // links to ancestor states
  void *links_des; // links to descendent states
} state;

typedef struct link
{
  DTYPE logp;
  state *ptr;
  struct link *next;
} link;

matrix* new_state(ITYPE n, DTYPE t, ITYPE *lins, DTYPE *values)
{
  state* stt;
  stt = (state *)malloc(sizeof(state));
  stt->n = n;
  stt->t = t;
  stt->lins = lins;
  stt->values = values;
  
  stt->links_anc = (link *)malloc(sizeof(link))
  stt->links_des = (link *)malloc(sizeof(link))
  return mat;
}

void print_state(state* stt)
{
  ITYPE n = stt->n;
  ITYPE t = stt->t;
  
  printf("state of %lu lineages at time %lf\n", n, t);
  
  ITYPE i, j;
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
      printf("%lf ", mat->data[i*n + j]);
    printf("\n");
  }
}



int main()
{
  //matrix* mat = newMatrix(10);
  //int idx[3] = {2,5,7}; 
  //addSquare(mat, idx, 3, 9);
  printf("%d", fib_(10));
  return 0;
}
