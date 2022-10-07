#define PY_SSIZE_T_CLEAN
#include <Python/Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

typedef double DTYPE;
typedef unsigned long ITYPE;

typedef struct list // list of pointers
{
  void **data;
  ITYPE used;
  ITYPE size;
} list;

list *new_list()
{
  list *lst = (list *)malloc(sizeof(list));
  lst->data = malloc(4 * sizeof(void *));
  lst->used = 0;
  lst->size = 4;
  return lst;
}

void insert_list(list *lst, void *ptr) 
{
  if (lst->used == lst->size)
  {
    lst->size *= 2;
    lst->data = realloc(lst->data, lst->size * sizeof(void *));
  }
  lst->data[lst->used++] = ptr;
}

void free_list(list *lst)
{
  free(lst->data);
  lst->data = NULL;
  lst->used = lst->size = 0;
}

/*
int main()
{
  list *lst = new_list();
  
  int i;
  void * ptr = NULL;
  for (i = 0; i < 100; i++)
    insert_list(lst, ptr);
  
  printf("%p\n", lst->data[9]);
  printf("%lu\n", lst->used);
  free_list(lst);
  
  return 0;
}
*/









typedef struct state
{
  ITYPE *values; // population labels of the lineages
  DTYPE logp;
  
  ITYPE num_anc; // number of ancestor states
  struct state **links_anc; // links to ancestor states
  DTYPE *logps_anc;
  
  ITYPE num_des; // number of descendent states
  struct state **links_des; // links to descendent states
  DTYPE *logps_des;
} state;

state* new_state(ITYPE *values)
{
  state* stt;
  stt = (state *)malloc(sizeof(state));
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




int main()
{
  ITYPE lins[5] = {4,9,2,10,31};
  ITYPE values[5] = {0,1,1,2,0};
  state *stt = new_state(5, 37.8, lins, values);
  print_state(stt);
  
  return 0;
}
