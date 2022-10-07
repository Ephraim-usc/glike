#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "numpy/arrayobject.h"

typedef double DTYPE;
typedef unsigned long ITYPE;

typedef struct state {
  ITYPE n; // number of lineages of the state
  DTYPE t; // time of the state
  ITYPE *lins; // lineages of the state
  DTYPE *values; // population labels of the lineages
  
  struct state **links_prev;
  struct state **links_next;
} state;

matrix* new_state(ITYPE n)
{
  state* stt;
  stt = (state *)malloc(sizeof(state));
  mat->n = n;
  mat->data = calloc(n*n,sizeof(DTYPE));
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
