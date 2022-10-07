#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

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



typedef struct hlist 
{
  struct hlist *next;
  char *key;
} hlist;

typedef struct htable
{
  ITYPE size;
  struct hlist **hlsts;
} htable;

def new_htable(ITYPE size)
{
  htable *htbl = (htable *)malloc(sizeof(htable));
  htbl->size = size;
  htbl->hlsts = (hlist **)malloc(size * sizeof(hlist *));
}


ITYPE hash(char *key, ITYPE size)
{
  ITYPE val;
  for (val = 0; *key; key++)
    val = 31 * val + *key;
  return val % size;
}


hlist *lookup_htable(htable *htbl, char *key)
{
  ITYPE val = hash(key, htbl->size);
  hlist *hlst;
  
  for (hlst = htable->hlsts[val]; hlst != NULL; hlst = hlst->next)
    if (strcmp(key, hlst->key) == 0)
      return hlst;
  
  return NULL;
}

void insert_htable(htable *htbl, char *key)
{
    hlist *hlst;
    ITYPE val;
    if ((np = lookup(name)) == NULL) /* not found */
    { 
        np = (struct nlist *) malloc(sizeof(*np));
        if (np == NULL || (np->name = strdup(name)) == NULL)
          return NULL;
        hashval = hash(name);
        np->next = hashtab[hashval];
        hashtab[hashval] = np;
    } 
  else /* already there */
    free((void *) np->defn); /*free previous defn */
  if ((np->defn = strdup(defn)) == NULL)
    return NULL;
  return np;
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
