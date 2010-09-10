/*-
 * Copyright 1997-2003 John-Mark Gurney.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *	$Id: fib.c,v 1.31 2003/01/14 10:11:30 jmg Exp $
 *
 */

/* Nils Homer - modified */

#include <stdlib.h>
#include <limits.h>
#include "fmap_alloc.h"
#include "fmap_fibheap.h"

#define swap(type, a, b)		\
  do {			\
      type c;		\
      c = a;		\
      a = b;		\
      b = c;		\
  } while (0)		\

#define INT_BITS        (sizeof(int) * 8)

#define fmap_fibheap_element_destroy(x) free((x))

static inline int
ceillog2(unsigned int a)
{
  int oa;
  int i;
  int b;

  oa = a;
  b = INT_BITS / 2;
  i = 0;
  while (b) {
      i = (i << 1);
      if(a >= (1 << b)) {
          a /= (1 << b);
          i = i | 1;
      } else
        a &= (1 << b) - 1;
      b /= 2;
  }
  if((1 << i) == oa)
    return i;
  else
    return i + 1;
}

/*
 * Private Heap Functions
 */
static void
fmap_fibheap_deleteel(fmap_fibheap_t *h, fmap_fibheap_element_t *x);
static void
fmap_fibheap_initheap(fmap_fibheap_t *new);
static void
fmap_fibheap_destroyheap(fmap_fibheap_t *h);

/*
 * begin of private element fuctions
 */
static fmap_fibheap_element_t *
fmap_fibheap_extractminel(fmap_fibheap_t *h);
static void
fmap_fibheap_insertrootlist(fmap_fibheap_t *h, fmap_fibheap_element_t *x);
static void
fmap_fibheap_removerootlist(fmap_fibheap_t *h, fmap_fibheap_element_t *x);
static void
fmap_fibheap_consolidate(fmap_fibheap_t *h);
static void
fmap_fibheap_heaplink(fmap_fibheap_t *h, fmap_fibheap_element_t *y, fmap_fibheap_element_t *x);
static void
fmap_fibheap_cut(fmap_fibheap_t *h, fmap_fibheap_element_t *x, fmap_fibheap_element_t *y);
static void
fmap_fibheap_cascading_cut(fmap_fibheap_t *h, fmap_fibheap_element_t *y);

/*
 * begining of handling elements of fmap_fibheap_t
 */
static fmap_fibheap_element_t *
fmap_fibheap_element_newelem();
static void
fmap_fibheap_element_initelem(fmap_fibheap_element_t *e);
static void
fmap_fibheap_element_insertafter(fmap_fibheap_element_t *a, fmap_fibheap_element_t *b);
static inline void
fmap_fibheap_element_insertbefore(fmap_fibheap_element_t *a, fmap_fibheap_element_t *b);
static fmap_fibheap_element_t *
fmap_fibheap_element_remove(fmap_fibheap_element_t *x);
static void
fmap_fibheap_checkcons(fmap_fibheap_t *h);
static int
fmap_fibheap_compare(fmap_fibheap_t *h, fmap_fibheap_element_t *a, fmap_fibheap_element_t *b);
static int
fmap_fibheap_comparedata(fmap_fibheap_t *h, int key, void *data, fmap_fibheap_element_t *b);
static void
fmap_fibheap_insertel(fmap_fibheap_t *h, fmap_fibheap_element_t *x);

/*
 * Private Heap Functions
 */
static void
fmap_fibheap_deleteel(fmap_fibheap_t *h, fmap_fibheap_element_t *x)
{
  void *data;
  int key;

  data = x->fmap_fibheap_element_data;
  key = x->fmap_fibheap_element_key;

  if(0 == h->fmap_fibheap_keys)
    fmap_fibheap_replacedata(h, x, h->fmap_fibheap_neginf);
  else
    fmap_fibheap_replacekey(h, x, INT_MIN);
  if(fmap_fibheap_extractminel(h) != x) {
      /*
       * XXX - This should never happen as fmap_fibheap_replace should set it
       * to min.
       */
      abort();
  }

  x->fmap_fibheap_element_data = data;
  x->fmap_fibheap_element_key = key;
}

static void
fmap_fibheap_initheap(fmap_fibheap_t *new)
{
  new->fmap_fibheap_cmp_fnct = NULL;
  new->fmap_fibheap_neginf = NULL;
  new->fmap_fibheap_n = 0;
  new->fmap_fibheap_Dl = -1;
  new->fmap_fibheap_cons = NULL;
  new->fmap_fibheap_min = NULL;
  new->fmap_fibheap_root = NULL;
  new->fmap_fibheap_keys = 0;
}

static void
fmap_fibheap_destroyheap(fmap_fibheap_t *h)
{
  h->fmap_fibheap_cmp_fnct = NULL;
  h->fmap_fibheap_neginf = NULL;
  if(h->fmap_fibheap_cons != NULL)
    free(h->fmap_fibheap_cons);
  h->fmap_fibheap_cons = NULL;
  free(h);
}

/*
 * begin of private element fuctions
 */
static fmap_fibheap_element_t *
fmap_fibheap_extractminel(fmap_fibheap_t *h)
{
  fmap_fibheap_element_t *ret;
  fmap_fibheap_element_t *x, *y, *orig;

  ret = h->fmap_fibheap_min;

  orig = NULL;
  /* put all the children on the root list */
  /* for true consistancy, we should use fmap_fibheap_element_remove */
  for(x = ret->fmap_fibheap_element_child; x != orig && x != NULL;) {
      if(orig == NULL)
        orig = x;
      y = x->fmap_fibheap_element_right;
      x->fmap_fibheap_element_p = NULL;
      fmap_fibheap_insertrootlist(h, x);
      x = y;
  }
  /* remove minimum from root list */
  fmap_fibheap_removerootlist(h, ret);
  h->fmap_fibheap_n--;

  /* if we aren't empty, consolidate the heap */
  if(h->fmap_fibheap_n == 0)
    h->fmap_fibheap_min = NULL;
  else {
      h->fmap_fibheap_min = ret->fmap_fibheap_element_right;
      fmap_fibheap_consolidate(h);
  }

  return ret;
}

static void
fmap_fibheap_insertrootlist(fmap_fibheap_t *h, fmap_fibheap_element_t *x)
{
  if(h->fmap_fibheap_root == NULL) {
      h->fmap_fibheap_root = x;
      x->fmap_fibheap_element_left = x;
      x->fmap_fibheap_element_right = x;
      return;
  }

  fmap_fibheap_element_insertafter(h->fmap_fibheap_root, x);
}

static void
fmap_fibheap_removerootlist(fmap_fibheap_t *h, fmap_fibheap_element_t *x)
{
  if(x->fmap_fibheap_element_left == x)
    h->fmap_fibheap_root = NULL;
  else
    h->fmap_fibheap_root = fmap_fibheap_element_remove(x);
}

static void
fmap_fibheap_consolidate(fmap_fibheap_t *h)
{
  fmap_fibheap_element_t **a;
  fmap_fibheap_element_t *w;
  fmap_fibheap_element_t *y;
  fmap_fibheap_element_t *x;
  int i;
  int d;
  int D;

  fmap_fibheap_checkcons(h);

  /* assign a the value of h->fmap_fibheap_cons so I don't have to rewrite code */
  D = h->fmap_fibheap_Dl + 1;
  a = h->fmap_fibheap_cons;

  for (i = 0; i < D; i++)
    a[i] = NULL;

  while ((w = h->fmap_fibheap_root) != NULL) {
      x = w;
      fmap_fibheap_removerootlist(h, w);
      d = x->fmap_fibheap_element_degree;
      /* XXX - assert that d < D */
      while(a[d] != NULL) {
          y = a[d];
          if(fmap_fibheap_compare(h, x, y) > 0)
            swap(fmap_fibheap_element_t *, x, y);
          fmap_fibheap_heaplink(h, y, x);
          a[d] = NULL;
          d++;
      }
      a[d] = x;
  }
  h->fmap_fibheap_min = NULL;
  for (i = 0; i < D; i++)
    if(a[i] != NULL) {
        fmap_fibheap_insertrootlist(h, a[i]);
        if(h->fmap_fibheap_min == NULL || fmap_fibheap_compare(h, a[i],
                                            h->fmap_fibheap_min) < 0)
          h->fmap_fibheap_min = a[i];
    }
}

static void
fmap_fibheap_heaplink(fmap_fibheap_t *h, fmap_fibheap_element_t *y, fmap_fibheap_element_t *x)
{
  /* make y a child of x */
  if(x->fmap_fibheap_element_child == NULL)
    x->fmap_fibheap_element_child = y;
  else
    fmap_fibheap_element_insertbefore(x->fmap_fibheap_element_child, y);
  y->fmap_fibheap_element_p = x;
  x->fmap_fibheap_element_degree++;
  y->fmap_fibheap_element_mark = 0;
}

static void
fmap_fibheap_cut(fmap_fibheap_t *h, fmap_fibheap_element_t *x, fmap_fibheap_element_t *y)
{
  fmap_fibheap_element_remove(x);
  y->fmap_fibheap_element_degree--;
  fmap_fibheap_insertrootlist(h, x);
  x->fmap_fibheap_element_p = NULL;
  x->fmap_fibheap_element_mark = 0;
}

static void
fmap_fibheap_cascading_cut(fmap_fibheap_t *h, fmap_fibheap_element_t *y)
{
  fmap_fibheap_element_t *z;

  while ((z = y->fmap_fibheap_element_p) != NULL) {
      if(y->fmap_fibheap_element_mark == 0) {
          y->fmap_fibheap_element_mark = 1;
          return;
      } else {
          fmap_fibheap_cut(h, y, z);
          y = z;
      }
  }
}

/*
 * begining of handling elements of fmap_fibheap_t
 */
static fmap_fibheap_element_t *
fmap_fibheap_element_newelem()
{
  fmap_fibheap_element_t *e;

  e = fmap_malloc(sizeof(fmap_fibheap_element_t), "e");
  fmap_fibheap_element_initelem(e);

  return e;
}

static void
fmap_fibheap_element_initelem(fmap_fibheap_element_t *e)
{
  e->fmap_fibheap_element_degree = 0;
  e->fmap_fibheap_element_mark = 0;
  e->fmap_fibheap_element_p = NULL;
  e->fmap_fibheap_element_child = NULL;
  e->fmap_fibheap_element_left = e;
  e->fmap_fibheap_element_right = e;
  e->fmap_fibheap_element_data = NULL;
}

static void
fmap_fibheap_element_insertafter(fmap_fibheap_element_t *a, fmap_fibheap_element_t *b)
{
  if(a == a->fmap_fibheap_element_right) {
      a->fmap_fibheap_element_right = b;
      a->fmap_fibheap_element_left = b;
      b->fmap_fibheap_element_right = a;
      b->fmap_fibheap_element_left = a;
  } else {
      b->fmap_fibheap_element_right = a->fmap_fibheap_element_right;
      a->fmap_fibheap_element_right->fmap_fibheap_element_left = b;
      a->fmap_fibheap_element_right = b;
      b->fmap_fibheap_element_left = a;
  }
}

static inline void
fmap_fibheap_element_insertbefore(fmap_fibheap_element_t *a, fmap_fibheap_element_t *b)
{
  fmap_fibheap_element_insertafter(a->fmap_fibheap_element_left, b);
}

static fmap_fibheap_element_t *
fmap_fibheap_element_remove(fmap_fibheap_element_t *x)
{
  fmap_fibheap_element_t *ret;

  if(x == x->fmap_fibheap_element_left)
    ret = NULL;
  else
    ret = x->fmap_fibheap_element_left;

  /* fix the parent pointer */
  if(x->fmap_fibheap_element_p != NULL && x->fmap_fibheap_element_p->fmap_fibheap_element_child == x)
    x->fmap_fibheap_element_p->fmap_fibheap_element_child = ret;

  x->fmap_fibheap_element_right->fmap_fibheap_element_left = x->fmap_fibheap_element_left;
  x->fmap_fibheap_element_left->fmap_fibheap_element_right = x->fmap_fibheap_element_right;

  /* clear out hanging pointers */
  x->fmap_fibheap_element_p = NULL;
  x->fmap_fibheap_element_left = x;
  x->fmap_fibheap_element_right = x;

  return ret;
}

static void
fmap_fibheap_checkcons(fmap_fibheap_t *h)
{
  int oDl;

  /* make sure we have enough memory allocated to "reorganize" */
  if(h->fmap_fibheap_Dl == -1 || h->fmap_fibheap_n > (1 << h->fmap_fibheap_Dl)) {
      oDl = h->fmap_fibheap_Dl;
      if((h->fmap_fibheap_Dl = ceillog2(h->fmap_fibheap_n) + 1) < 8)
        h->fmap_fibheap_Dl = 8;
      if(oDl != h->fmap_fibheap_Dl)
        h->fmap_fibheap_cons = fmap_realloc(h->fmap_fibheap_cons, sizeof(fmap_fibheap_element_t) * (h->fmap_fibheap_Dl + 1), "h->fmap_fibheap_cons");
      if(h->fmap_fibheap_cons == NULL)
        abort();
  }
}

static int
fmap_fibheap_compare(fmap_fibheap_t *h, fmap_fibheap_element_t *a, fmap_fibheap_element_t *b)
{
  if(1 == h->fmap_fibheap_keys) {
      if(a->fmap_fibheap_element_key < b->fmap_fibheap_element_key)
        return -1;
      if(a->fmap_fibheap_element_key == b->fmap_fibheap_element_key)
        return 0;
      return 1;
  } else {
    return h->fmap_fibheap_cmp_fnct(a->fmap_fibheap_element_data, b->fmap_fibheap_element_data);
  }
}

static int
fmap_fibheap_comparedata(fmap_fibheap_t *h, int key, void *data, fmap_fibheap_element_t *b)
{
  fmap_fibheap_element_t a;

  a.fmap_fibheap_element_key = key;
  a.fmap_fibheap_element_data = data;

  return fmap_fibheap_compare(h, &a, b);
}

static void
fmap_fibheap_insertel(fmap_fibheap_t *h, fmap_fibheap_element_t *x)
{
  fmap_fibheap_insertrootlist(h, x);

  if(h->fmap_fibheap_min == NULL 
     || (1 == h->fmap_fibheap_keys 
         && x->fmap_fibheap_element_key < h->fmap_fibheap_min->fmap_fibheap_element_key) 
     || (0 == h->fmap_fibheap_keys 
         && h->fmap_fibheap_cmp_fnct(x->fmap_fibheap_element_data, 
                                     h->fmap_fibheap_min->fmap_fibheap_element_data) < 0)) {
    h->fmap_fibheap_min = x;
  }

  h->fmap_fibheap_n++;

}

/*
 * Public Heap Functions
 */
fmap_fibheap_t *
fmap_fibheap_makekeyheap()
{
  fmap_fibheap_t *n;

  n = fmap_malloc(sizeof(fmap_fibheap_t), "n");

  fmap_fibheap_initheap(n);
  n->fmap_fibheap_keys = 1;

  return n;
}

fmap_fibheap_t *
fmap_fibheap_makeheap(fmap_fibheap_voidcmp fnct)
{
  fmap_fibheap_t *n;

  n = fmap_malloc(sizeof(fmap_fibheap_t), "n");

  fmap_fibheap_initheap(n);
  fmap_fibheap_setcmp(n, fnct);
  n->fmap_fibheap_keys = 0;

  return n;
}

fmap_fibheap_voidcmp
fmap_fibheap_setcmp(fmap_fibheap_t *h, fmap_fibheap_voidcmp fnct)
{
  fmap_fibheap_voidcmp oldfnct;

  if(1 == h->fmap_fibheap_keys) fmap_error("void heap required", Exit, OutOfRange);
  oldfnct = h->fmap_fibheap_cmp_fnct;
  h->fmap_fibheap_cmp_fnct = fnct;

  return oldfnct;
}

void *
fmap_fibheap_setneginf(fmap_fibheap_t *h, void *data)
{
  void *old;

  if(1 == h->fmap_fibheap_keys) fmap_error("void heap required", Exit, OutOfRange);
  old = h->fmap_fibheap_neginf;
  h->fmap_fibheap_neginf = data;

  return old;
}

fmap_fibheap_t *
fmap_fibheap_union(fmap_fibheap_t *ha, fmap_fibheap_t *hb)
{
  fmap_fibheap_element_t *x;

  if(ha->fmap_fibheap_root == NULL || hb->fmap_fibheap_root == NULL) {
      /* either one or both are empty */
      if(ha->fmap_fibheap_root == NULL) {
          fmap_fibheap_destroyheap(ha);
          return hb;
      } else {
          fmap_fibheap_destroyheap(hb);
          return ha;
      }
  }
  ha->fmap_fibheap_root->fmap_fibheap_element_left->fmap_fibheap_element_right = hb->fmap_fibheap_root;
  hb->fmap_fibheap_root->fmap_fibheap_element_left->fmap_fibheap_element_right = ha->fmap_fibheap_root;
  x = ha->fmap_fibheap_root->fmap_fibheap_element_left;
  ha->fmap_fibheap_root->fmap_fibheap_element_left = hb->fmap_fibheap_root->fmap_fibheap_element_left;
  hb->fmap_fibheap_root->fmap_fibheap_element_left = x;
  ha->fmap_fibheap_n += hb->fmap_fibheap_n;
  /*
   * we probably should also keep stats on number of unions
   */

  /* set fmap_fibheap_min if necessary */
  if(fmap_fibheap_compare(ha, hb->fmap_fibheap_min, ha->fmap_fibheap_min) < 0)
    ha->fmap_fibheap_min = hb->fmap_fibheap_min;

  fmap_fibheap_destroyheap(hb);
  return ha;
}

void
fmap_fibheap_deleteheap(fmap_fibheap_t *h)
{
  /*
   * We could do this even faster by walking each binomial tree, but
   * this is simpler to code.
   */
  while (h->fmap_fibheap_min != NULL)
    fmap_fibheap_element_destroy(fmap_fibheap_extractminel(h));

  fmap_fibheap_destroyheap(h);
}

/*
 * Public Key Heap Functions
 */
fmap_fibheap_element_t *
fmap_fibheap_insertkey(fmap_fibheap_t *h, int key, void *data)
{
  fmap_fibheap_element_t *x;

  if(0 == h->fmap_fibheap_keys) fmap_error("key heap required", Exit, OutOfRange);

  if((x = fmap_fibheap_element_newelem()) == NULL)
    return NULL;

  /* just insert on root list, and make sure it's not the new min */
  x->fmap_fibheap_element_data = data;
  x->fmap_fibheap_element_key = key;

  fmap_fibheap_insertel(h, x);

  return x;
}

int
fmap_fibheap_minkey(fmap_fibheap_t *h)
{
  if(0 == h->fmap_fibheap_keys) fmap_error("key heap required", Exit, OutOfRange);
  if(h->fmap_fibheap_min == NULL)
    return INT_MIN;
  return h->fmap_fibheap_min->fmap_fibheap_element_key;
}

int
fmap_fibheap_replacekey(fmap_fibheap_t *h, fmap_fibheap_element_t *x, int key)
{
  int ret;

  if(0 == h->fmap_fibheap_keys) fmap_error("key heap required", Exit, OutOfRange);
  ret = x->fmap_fibheap_element_key;
  (void)fmap_fibheap_replacekeydata(h, x, key, x->fmap_fibheap_element_data);

  return ret;
}

void *
fmap_fibheap_replacekeydata(fmap_fibheap_t *h, fmap_fibheap_element_t *x, int key, void *data)
{
  void *odata;
  int okey;
  fmap_fibheap_element_t *y;
  int r;

  if(0 == h->fmap_fibheap_keys) fmap_error("key heap required", Exit, OutOfRange);

  odata = x->fmap_fibheap_element_data;
  okey = x->fmap_fibheap_element_key;

  /*
   * we can increase a key by deleting and reinserting, that
   * requires O(lgn) time.
   */
  if((r = fmap_fibheap_comparedata(h, key, data, x)) > 0) {
      /* XXX - bad code! */
      abort();
      fmap_fibheap_deleteel(h, x);

      x->fmap_fibheap_element_data = data;
      x->fmap_fibheap_element_key = key;

      fmap_fibheap_insertel(h, x);

      return odata;
  }

  x->fmap_fibheap_element_data = data;
  x->fmap_fibheap_element_key = key;

  /* because they are equal, we don't have to do anything */
  if(r == 0)
    return odata;

  y = x->fmap_fibheap_element_p;

  if(1 == h->fmap_fibheap_keys && okey == key)
    return odata;

  if(y != NULL && fmap_fibheap_compare(h, x, y) <= 0) {
      fmap_fibheap_cut(h, x, y);
      fmap_fibheap_cascading_cut(h, y);
  }

  /*
   * the = is so that the call from fmap_fibheap_delete will delete the proper
   * element.
   */
  if(fmap_fibheap_compare(h, x, h->fmap_fibheap_min) <= 0)
    h->fmap_fibheap_min = x;

  return odata;
}

/*
 * Public void * Heap Functions
 */
/*
 * this will return these values:
 *	NULL	failed for some reason
 *	ptr	token to use for manipulation of data
 */
fmap_fibheap_element_t *
fmap_fibheap_insert(fmap_fibheap_t *h, void *data)
{
  fmap_fibheap_element_t *x;

  if(1 == h->fmap_fibheap_keys) fmap_error("void heap required", Exit, OutOfRange);
  if((x = fmap_fibheap_element_newelem()) == NULL)
    return NULL;

  /* just insert on root list, and make sure it's not the new min */
  x->fmap_fibheap_element_data = data;

  fmap_fibheap_insertel(h, x);

  return x;
}

void *
fmap_fibheap_min(fmap_fibheap_t *h)
{
  if(h->fmap_fibheap_min == NULL)
    return NULL;
  return h->fmap_fibheap_min->fmap_fibheap_element_data;
}

void *
fmap_fibheap_extractmin(fmap_fibheap_t *h)
{
  fmap_fibheap_element_t *z;
  void *ret;

  ret = NULL;

  if(h->fmap_fibheap_min != NULL) {
      z = fmap_fibheap_extractminel(h);
      ret = z->fmap_fibheap_element_data;
#ifndef FMAP_FIBHEAP_NO_FREE
      fmap_fibheap_element_destroy(z);
#endif

  }

  return ret;
}

void *
fmap_fibheap_replacedata(fmap_fibheap_t *h, fmap_fibheap_element_t *x, void *data)
{
  return fmap_fibheap_replacekeydata(h, x, x->fmap_fibheap_element_key, data);
}

void *
fmap_fibheap_delete(fmap_fibheap_t *h, fmap_fibheap_element_t *x)
{
  void *k;

  k = x->fmap_fibheap_element_data;
  if(0 == h->fmap_fibheap_keys)
    fmap_fibheap_replacedata(h, x, h->fmap_fibheap_neginf);
  else
    fmap_fibheap_replacekey(h, x, INT_MIN);
  fmap_fibheap_extractmin(h);

  return k;
}

