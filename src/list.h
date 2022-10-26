#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>



typedef struct hnode
{
  struct hnode *next;
  int *key;
  void *pointer;
} hnode;

typedef struct htable
{
  int size;
  hnode **hnodes;
} htable;

htable *new_htable(int size)
{
  htable *htbl = (htable *)malloc(sizeof(htable));
  htbl->size = size;
  htbl->hnodes = (hnode **)malloc(size * sizeof(hnode *));
  return htbl;
}

int hash(int *key, int len, int size)
{
  int *p = key;
  int val = 0;
  for (; p < key + len; p++)
    val = (31 * val + *p) % size;
  return val;
}

int compare(int *a, int *b, int len)
{
  int i = 0;
  while (*(a+i) == *(b+i) && i < len-1)
    i++;
  return (*(a+i) > *(b+i)) - (*(a+i) < *(b+i));
}

hnode *lookup_htable(htable *htbl, int *key, int len)
{
  int val = hash(key, len, htbl->size);
  hnode *hnd;
  
  for (hnd = htbl->hnodes[val]; hnd != NULL; hnd = hnd->next)
    if (compare(key, hnd->key, len) == 0)
      return hnd;
  
  return NULL;
}

hnode *insert_htable(htable *htbl, int *key, int len)
{
  int val = hash(key, len, htbl->size);
  hnode *hnd;
  
  for (hnd = htbl->hnodes[val]; hnd != NULL; hnd = hnd->next)
    if (compare(key, hnd->key, len) == 0)
      return hnd;
  
  hnd = (hnode *) malloc(sizeof(hnode));
  hnd->key = key;
  hnd->pointer = NULL;
  
  hnd->next = htbl->hnodes[val];
  htbl->hnodes[val] = hnd;
  return hnd;
}

/*
int main()
{
  int a[5] = {0,5,2,1,4};
  int b[5] = {0,5,4,1,4};
  int c[5] = {0,6,2,1,4};
  int d[5] = {0,5,2,1,35};
  int e[5] = {0,5,4,1,4};

  htable *htbl = new_htable(10);
  printf("%d\n", hash(a, 5, 10));
  printf("%p\n", insert_htable(htbl, a, 5));
  printf("%d\n", hash(b, 5, 10));
  printf("%p\n", insert_htable(htbl, b, 5));
  printf("%d\n", hash(c, 5, 10));
  printf("%p\n", insert_htable(htbl, c, 5));
  printf("%d\n", hash(d, 5, 10));
  printf("%p\n", insert_htable(htbl, d, 5));
  printf("%d\n", hash(e, 5, 10));
  printf("%p\n", insert_htable(htbl, e, 5));
  
  return 0;
}
*/


























typedef double DTYPE;
typedef unsigned long ITYPE;


typedef struct list // list of DTYPE numbers
{
  DTYPE *data;
  ITYPE used;
  ITYPE size;
} list;

list *new_list()
{
  list *lst = (list *)malloc(sizeof(list));
  lst->data = malloc(1 * sizeof(DTYPE));
  lst->used = 0;
  lst->size = 1;
  return lst;
}

void insert_list(list *lst, DTYPE a) 
{
  if (lst->used == lst->size)
  {
    lst->size *= 2;
    lst->data = realloc(lst->data, lst->size * sizeof(DTYPE));
  }
  lst->data[lst->used++] = a;
}

void free_list(list *lst)
{
  free(lst->data);
  lst->data = NULL;
  lst->used = lst->size = 0;
}


typedef struct plist // list of pointers
{
  void **data;
  ITYPE used;
  ITYPE size;
} plist;

plist *new_plist()
{
  plist *lst = (plist *)malloc(sizeof(plist));
  lst->data = malloc(4 * sizeof(void *));
  lst->used = 0;
  lst->size = 4;
  return lst;
}

void insert_plist(plist *lst, void *ptr) 
{
  if (lst->used == lst->size)
  {
    lst->size *= 2;
    lst->data = realloc(lst->data, lst->size * sizeof(void *));
  }
  lst->data[lst->used++] = ptr;
}

void free_plist(plist *lst)
{
  free(lst->data);
  lst->data = NULL;
  lst->used = lst->size = 0;
}

/*
int main()
{
  DTYPE a = 0.36;
  void * ptr = NULL;
  int i;
  
  
  list *lst = new_list();
  for (i = 0; i < 100; i++)
    insert_list(lst, a);
  
  printf("%lf\n", lst->data[9]);
  printf("%lu\n", lst->used);
  free_list(lst);
  
  
  plist *plst = new_plist();
  for (i = 0; i < 100; i++)
    insert_plist(plst, ptr);
  
  printf("%p\n", plst->data[9]);
  printf("%lu\n", plst->used);
  free_plist(plst);
  
  
  return 0;
}
*/














