#define PY_SSIZE_T_CLEAN
#include <Python/Python.h>
//#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

typedef double DTYPE;
typedef unsigned long ITYPE;

typedef struct state
{
  ITYPE n; // number of lineages of the state
  DTYPE t; // time of the state
  ITYPE *lins; // lineages of the state
  ITYPE *values; // population labels of the lineages
  DTYPE logp;
  
  ITYPE num_anc; // number of ancestor states
  void *links_anc; // links to ancestor states
  DTYPE *logps_anc;
  
  ITYPE num_des; // number of descendent states
  void *links_des; // links to descendent states
  DTYPE *logps_des;
} state;

state* new_state(ITYPE n, DTYPE t, ITYPE *lins, ITYPE *values)
{
  state* stt;
  stt = (state *)malloc(sizeof(state));
  stt->n = n;
  stt->t = t;
  stt->lins = lins;
  stt->values = values;
  
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
  
  for(i = 0; i < stt->num_anc; i++)
    printf("%lf", *(stt->logps_anc+i));
  printf("\n");
  
  for(i = 0; i < stt->num_des; i++)
    printf("%lf", *(stt->logps_des+i));
  printf("\n");
}



int main()
{
  ITYPE lins[5] = {4,9,2,10,31};
  ITYPE values[5] = {0,1,1,2,0};
  state *stt = new_state(5, 37.8, lins, values);
  print_state(stt);
  
  return 0;
}
