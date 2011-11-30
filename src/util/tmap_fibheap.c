/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
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
#include "tmap_alloc.h"
#include "tmap_fibheap.h"

#define swap(type, a, b)		\
  do {			\
      type c;		\
      c = a;		\
      a = b;		\
      b = c;		\
  } while (0)		\

#define INT_BITS        (sizeof(int) * 8)

#define tmap_fibheap_element_destroy(x) free((x))

static inline int
ceillog2(int a)
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
tmap_fibheap_deleteel(tmap_fibheap_t *h, tmap_fibheap_element_t *x);
static void
tmap_fibheap_initheap(tmap_fibheap_t *new);
static void
tmap_fibheap_destroyheap(tmap_fibheap_t *h);

/*
 * begin of private element fuctions
 */
static tmap_fibheap_element_t *
tmap_fibheap_extractminel(tmap_fibheap_t *h);
static void
tmap_fibheap_insertrootlist(tmap_fibheap_t *h, tmap_fibheap_element_t *x);
static void
tmap_fibheap_removerootlist(tmap_fibheap_t *h, tmap_fibheap_element_t *x);
static void
tmap_fibheap_consolidate(tmap_fibheap_t *h);
static void
tmap_fibheap_heaplink(tmap_fibheap_element_t *y, tmap_fibheap_element_t *x);
static void
tmap_fibheap_cut(tmap_fibheap_t *h, tmap_fibheap_element_t *x, tmap_fibheap_element_t *y);
static void
tmap_fibheap_cascading_cut(tmap_fibheap_t *h, tmap_fibheap_element_t *y);

/*
 * begining of handling elements of tmap_fibheap_t
 */
static tmap_fibheap_element_t *
tmap_fibheap_element_newelem();
static void
tmap_fibheap_element_initelem(tmap_fibheap_element_t *e);
static void
tmap_fibheap_element_insertafter(tmap_fibheap_element_t *a, tmap_fibheap_element_t *b);
static inline void
tmap_fibheap_element_insertbefore(tmap_fibheap_element_t *a, tmap_fibheap_element_t *b);
static tmap_fibheap_element_t *
tmap_fibheap_element_remove(tmap_fibheap_element_t *x);
static void
tmap_fibheap_checkcons(tmap_fibheap_t *h);
static int
tmap_fibheap_compare(tmap_fibheap_t *h, tmap_fibheap_element_t *a, tmap_fibheap_element_t *b);
static int
tmap_fibheap_comparedata(tmap_fibheap_t *h, int key, void *data, tmap_fibheap_element_t *b);
static void
tmap_fibheap_insertel(tmap_fibheap_t *h, tmap_fibheap_element_t *x);

/*
 * Private Heap Functions
 */
static void
tmap_fibheap_deleteel(tmap_fibheap_t *h, tmap_fibheap_element_t *x)
{
  void *data;
  int key;

  data = x->tmap_fibheap_element_data;
  key = x->tmap_fibheap_element_key;

  if(0 == h->tmap_fibheap_keys)
    tmap_fibheap_replacedata(h, x, h->tmap_fibheap_neginf);
  else
    tmap_fibheap_replacekey(h, x, INT_MIN);
  if(tmap_fibheap_extractminel(h) != x) {
      /*
       * XXX - This should never happen as tmap_fibheap_replace should set it
       * to min.
       */
      abort();
  }

  x->tmap_fibheap_element_data = data;
  x->tmap_fibheap_element_key = key;
}

static void
tmap_fibheap_initheap(tmap_fibheap_t *new)
{
  new->tmap_fibheap_cmp_fnct = NULL;
  new->tmap_fibheap_neginf = NULL;
  new->tmap_fibheap_n = 0;
  new->tmap_fibheap_Dl = -1;
  new->tmap_fibheap_cons = NULL;
  new->tmap_fibheap_min = NULL;
  new->tmap_fibheap_root = NULL;
  new->tmap_fibheap_keys = 0;
}

static void
tmap_fibheap_destroyheap(tmap_fibheap_t *h)
{
  h->tmap_fibheap_cmp_fnct = NULL;
  h->tmap_fibheap_neginf = NULL;
  if(h->tmap_fibheap_cons != NULL)
    free(h->tmap_fibheap_cons);
  h->tmap_fibheap_cons = NULL;
  free(h);
}

/*
 * begin of private element fuctions
 */
static tmap_fibheap_element_t *
tmap_fibheap_extractminel(tmap_fibheap_t *h)
{
  tmap_fibheap_element_t *ret;
  tmap_fibheap_element_t *x, *y, *orig;

  ret = h->tmap_fibheap_min;

  orig = NULL;
  /* put all the children on the root list */
  /* for true consistancy, we should use tmap_fibheap_element_remove */
  for(x = ret->tmap_fibheap_element_child; x != orig && x != NULL;) {
      if(orig == NULL)
        orig = x;
      y = x->tmap_fibheap_element_right;
      x->tmap_fibheap_element_p = NULL;
      tmap_fibheap_insertrootlist(h, x);
      x = y;
  }
  /* remove minimum from root list */
  tmap_fibheap_removerootlist(h, ret);
  h->tmap_fibheap_n--;

  /* if we aren't empty, consolidate the heap */
  if(h->tmap_fibheap_n == 0)
    h->tmap_fibheap_min = NULL;
  else {
      h->tmap_fibheap_min = ret->tmap_fibheap_element_right;
      tmap_fibheap_consolidate(h);
  }

  return ret;
}

static void
tmap_fibheap_insertrootlist(tmap_fibheap_t *h, tmap_fibheap_element_t *x)
{
  if(h->tmap_fibheap_root == NULL) {
      h->tmap_fibheap_root = x;
      x->tmap_fibheap_element_left = x;
      x->tmap_fibheap_element_right = x;
      return;
  }

  tmap_fibheap_element_insertafter(h->tmap_fibheap_root, x);
}

static void
tmap_fibheap_removerootlist(tmap_fibheap_t *h, tmap_fibheap_element_t *x)
{
  if(x->tmap_fibheap_element_left == x)
    h->tmap_fibheap_root = NULL;
  else
    h->tmap_fibheap_root = tmap_fibheap_element_remove(x);
}

static void
tmap_fibheap_consolidate(tmap_fibheap_t *h)
{
  tmap_fibheap_element_t **a;
  tmap_fibheap_element_t *w;
  tmap_fibheap_element_t *y;
  tmap_fibheap_element_t *x;
  int i;
  int d;
  int D;

  tmap_fibheap_checkcons(h);

  /* assign a the value of h->tmap_fibheap_cons so I don't have to rewrite code */
  D = h->tmap_fibheap_Dl + 1;
  a = h->tmap_fibheap_cons;

  for (i = 0; i < D; i++)
    a[i] = NULL;

  while ((w = h->tmap_fibheap_root) != NULL) {
      x = w;
      tmap_fibheap_removerootlist(h, w);
      d = x->tmap_fibheap_element_degree;
      /* XXX - assert that d < D */
      while(a[d] != NULL) {
          y = a[d];
          if(tmap_fibheap_compare(h, x, y) > 0)
            swap(tmap_fibheap_element_t *, x, y);
          tmap_fibheap_heaplink(y, x);
          a[d] = NULL;
          d++;
      }
      a[d] = x;
  }
  h->tmap_fibheap_min = NULL;
  for (i = 0; i < D; i++)
    if(a[i] != NULL) {
        tmap_fibheap_insertrootlist(h, a[i]);
        if(h->tmap_fibheap_min == NULL || tmap_fibheap_compare(h, a[i],
                                            h->tmap_fibheap_min) < 0)
          h->tmap_fibheap_min = a[i];
    }
}

static void
tmap_fibheap_heaplink(tmap_fibheap_element_t *y, tmap_fibheap_element_t *x)
{
  /* make y a child of x */
  if(x->tmap_fibheap_element_child == NULL)
    x->tmap_fibheap_element_child = y;
  else
    tmap_fibheap_element_insertbefore(x->tmap_fibheap_element_child, y);
  y->tmap_fibheap_element_p = x;
  x->tmap_fibheap_element_degree++;
  y->tmap_fibheap_element_mark = 0;
}

static void
tmap_fibheap_cut(tmap_fibheap_t *h, tmap_fibheap_element_t *x, tmap_fibheap_element_t *y)
{
  tmap_fibheap_element_remove(x);
  y->tmap_fibheap_element_degree--;
  tmap_fibheap_insertrootlist(h, x);
  x->tmap_fibheap_element_p = NULL;
  x->tmap_fibheap_element_mark = 0;
}

static void
tmap_fibheap_cascading_cut(tmap_fibheap_t *h, tmap_fibheap_element_t *y)
{
  tmap_fibheap_element_t *z;

  while ((z = y->tmap_fibheap_element_p) != NULL) {
      if(y->tmap_fibheap_element_mark == 0) {
          y->tmap_fibheap_element_mark = 1;
          return;
      } else {
          tmap_fibheap_cut(h, y, z);
          y = z;
      }
  }
}

/*
 * begining of handling elements of tmap_fibheap_t
 */
static tmap_fibheap_element_t *
tmap_fibheap_element_newelem()
{
  tmap_fibheap_element_t *e;

  e = tmap_malloc(sizeof(tmap_fibheap_element_t), "e");
  tmap_fibheap_element_initelem(e);

  return e;
}

static void
tmap_fibheap_element_initelem(tmap_fibheap_element_t *e)
{
  e->tmap_fibheap_element_degree = 0;
  e->tmap_fibheap_element_mark = 0;
  e->tmap_fibheap_element_p = NULL;
  e->tmap_fibheap_element_child = NULL;
  e->tmap_fibheap_element_left = e;
  e->tmap_fibheap_element_right = e;
  e->tmap_fibheap_element_data = NULL;
}

static void
tmap_fibheap_element_insertafter(tmap_fibheap_element_t *a, tmap_fibheap_element_t *b)
{
  if(a == a->tmap_fibheap_element_right) {
      a->tmap_fibheap_element_right = b;
      a->tmap_fibheap_element_left = b;
      b->tmap_fibheap_element_right = a;
      b->tmap_fibheap_element_left = a;
  } else {
      b->tmap_fibheap_element_right = a->tmap_fibheap_element_right;
      a->tmap_fibheap_element_right->tmap_fibheap_element_left = b;
      a->tmap_fibheap_element_right = b;
      b->tmap_fibheap_element_left = a;
  }
}

static inline void
tmap_fibheap_element_insertbefore(tmap_fibheap_element_t *a, tmap_fibheap_element_t *b)
{
  tmap_fibheap_element_insertafter(a->tmap_fibheap_element_left, b);
}

static tmap_fibheap_element_t *
tmap_fibheap_element_remove(tmap_fibheap_element_t *x)
{
  tmap_fibheap_element_t *ret;

  if(x == x->tmap_fibheap_element_left)
    ret = NULL;
  else
    ret = x->tmap_fibheap_element_left;

  /* fix the parent pointer */
  if(x->tmap_fibheap_element_p != NULL && x->tmap_fibheap_element_p->tmap_fibheap_element_child == x)
    x->tmap_fibheap_element_p->tmap_fibheap_element_child = ret;

  x->tmap_fibheap_element_right->tmap_fibheap_element_left = x->tmap_fibheap_element_left;
  x->tmap_fibheap_element_left->tmap_fibheap_element_right = x->tmap_fibheap_element_right;

  /* clear out hanging pointers */
  x->tmap_fibheap_element_p = NULL;
  x->tmap_fibheap_element_left = x;
  x->tmap_fibheap_element_right = x;

  return ret;
}

static void
tmap_fibheap_checkcons(tmap_fibheap_t *h)
{
  int oDl;

  /* make sure we have enough memory allocated to "reorganize" */
  if(h->tmap_fibheap_Dl == -1 || h->tmap_fibheap_n > (1 << h->tmap_fibheap_Dl)) {
      oDl = h->tmap_fibheap_Dl;
      if((h->tmap_fibheap_Dl = ceillog2(h->tmap_fibheap_n) + 1) < 8)
        h->tmap_fibheap_Dl = 8;
      if(oDl != h->tmap_fibheap_Dl)
        h->tmap_fibheap_cons = tmap_realloc(h->tmap_fibheap_cons, sizeof(tmap_fibheap_element_t) * (h->tmap_fibheap_Dl + 1), "h->tmap_fibheap_cons");
      if(h->tmap_fibheap_cons == NULL)
        abort();
  }
}

static int
tmap_fibheap_compare(tmap_fibheap_t *h, tmap_fibheap_element_t *a, tmap_fibheap_element_t *b)
{
  if(1 == h->tmap_fibheap_keys) {
      if(a->tmap_fibheap_element_key < b->tmap_fibheap_element_key)
        return -1;
      if(a->tmap_fibheap_element_key == b->tmap_fibheap_element_key)
        return 0;
      return 1;
  } else {
    return h->tmap_fibheap_cmp_fnct(a->tmap_fibheap_element_data, b->tmap_fibheap_element_data);
  }
}

static int
tmap_fibheap_comparedata(tmap_fibheap_t *h, int key, void *data, tmap_fibheap_element_t *b)
{
  tmap_fibheap_element_t a;

  a.tmap_fibheap_element_key = key;
  a.tmap_fibheap_element_data = data;

  return tmap_fibheap_compare(h, &a, b);
}

static void
tmap_fibheap_insertel(tmap_fibheap_t *h, tmap_fibheap_element_t *x)
{
  tmap_fibheap_insertrootlist(h, x);

  if(h->tmap_fibheap_min == NULL 
     || (1 == h->tmap_fibheap_keys 
         && x->tmap_fibheap_element_key < h->tmap_fibheap_min->tmap_fibheap_element_key) 
     || (0 == h->tmap_fibheap_keys 
         && h->tmap_fibheap_cmp_fnct(x->tmap_fibheap_element_data, 
                                     h->tmap_fibheap_min->tmap_fibheap_element_data) < 0)) {
    h->tmap_fibheap_min = x;
  }

  h->tmap_fibheap_n++;

}

/*
 * Public Heap Functions
 */
tmap_fibheap_t *
tmap_fibheap_makekeyheap()
{
  tmap_fibheap_t *n;

  n = tmap_malloc(sizeof(tmap_fibheap_t), "n");

  tmap_fibheap_initheap(n);
  n->tmap_fibheap_keys = 1;

  return n;
}

tmap_fibheap_t *
tmap_fibheap_makeheap(tmap_fibheap_voidcmp fnct)
{
  tmap_fibheap_t *n;

  n = tmap_malloc(sizeof(tmap_fibheap_t), "n");

  tmap_fibheap_initheap(n);
  tmap_fibheap_setcmp(n, fnct);
  n->tmap_fibheap_keys = 0;

  return n;
}

tmap_fibheap_voidcmp
tmap_fibheap_setcmp(tmap_fibheap_t *h, tmap_fibheap_voidcmp fnct)
{
  tmap_fibheap_voidcmp oldfnct;

  if(1 == h->tmap_fibheap_keys) tmap_error("void heap required", Exit, OutOfRange);
  oldfnct = h->tmap_fibheap_cmp_fnct;
  h->tmap_fibheap_cmp_fnct = fnct;

  return oldfnct;
}

void *
tmap_fibheap_setneginf(tmap_fibheap_t *h, void *data)
{
  void *old;

  if(1 == h->tmap_fibheap_keys) tmap_error("void heap required", Exit, OutOfRange);
  old = h->tmap_fibheap_neginf;
  h->tmap_fibheap_neginf = data;

  return old;
}

tmap_fibheap_t *
tmap_fibheap_union(tmap_fibheap_t *ha, tmap_fibheap_t *hb)
{
  tmap_fibheap_element_t *x;

  if(ha->tmap_fibheap_root == NULL || hb->tmap_fibheap_root == NULL) {
      /* either one or both are empty */
      if(ha->tmap_fibheap_root == NULL) {
          tmap_fibheap_destroyheap(ha);
          return hb;
      } else {
          tmap_fibheap_destroyheap(hb);
          return ha;
      }
  }
  ha->tmap_fibheap_root->tmap_fibheap_element_left->tmap_fibheap_element_right = hb->tmap_fibheap_root;
  hb->tmap_fibheap_root->tmap_fibheap_element_left->tmap_fibheap_element_right = ha->tmap_fibheap_root;
  x = ha->tmap_fibheap_root->tmap_fibheap_element_left;
  ha->tmap_fibheap_root->tmap_fibheap_element_left = hb->tmap_fibheap_root->tmap_fibheap_element_left;
  hb->tmap_fibheap_root->tmap_fibheap_element_left = x;
  ha->tmap_fibheap_n += hb->tmap_fibheap_n;
  /*
   * we probably should also keep stats on number of unions
   */

  /* set tmap_fibheap_min if necessary */
  if(tmap_fibheap_compare(ha, hb->tmap_fibheap_min, ha->tmap_fibheap_min) < 0)
    ha->tmap_fibheap_min = hb->tmap_fibheap_min;

  tmap_fibheap_destroyheap(hb);
  return ha;
}

void
tmap_fibheap_deleteheap(tmap_fibheap_t *h)
{
  /*
   * We could do this even faster by walking each binomial tree, but
   * this is simpler to code.
   */
  while (h->tmap_fibheap_min != NULL)
    tmap_fibheap_element_destroy(tmap_fibheap_extractminel(h));

  tmap_fibheap_destroyheap(h);
}

/*
 * Public Key Heap Functions
 */
tmap_fibheap_element_t *
tmap_fibheap_insertkey(tmap_fibheap_t *h, int key, void *data)
{
  tmap_fibheap_element_t *x;

  if(0 == h->tmap_fibheap_keys) tmap_error("key heap required", Exit, OutOfRange);

  if((x = tmap_fibheap_element_newelem()) == NULL) {
    return NULL;
  }

  /* just insert on root list, and make sure it's not the new min */
  x->tmap_fibheap_element_data = data;
  x->tmap_fibheap_element_key = key;

  tmap_fibheap_insertel(h, x);

  return x;
}

int
tmap_fibheap_minkey(tmap_fibheap_t *h)
{
  if(0 == h->tmap_fibheap_keys) tmap_error("key heap required", Exit, OutOfRange);
  if(h->tmap_fibheap_min == NULL)
    return INT_MIN;
  return h->tmap_fibheap_min->tmap_fibheap_element_key;
}

int
tmap_fibheap_replacekey(tmap_fibheap_t *h, tmap_fibheap_element_t *x, int key)
{
  int ret;

  if(0 == h->tmap_fibheap_keys) tmap_error("key heap required", Exit, OutOfRange);
  ret = x->tmap_fibheap_element_key;
  (void)tmap_fibheap_replacekeydata(h, x, key, x->tmap_fibheap_element_data);

  return ret;
}

void *
tmap_fibheap_replacekeydata(tmap_fibheap_t *h, tmap_fibheap_element_t *x, int key, void *data)
{
  void *odata;
  int okey;
  tmap_fibheap_element_t *y;
  int r;

  if(0 == h->tmap_fibheap_keys) tmap_error("key heap required", Exit, OutOfRange);

  odata = x->tmap_fibheap_element_data;
  okey = x->tmap_fibheap_element_key;

  /*
   * we can increase a key by deleting and reinserting, that
   * requires O(lgn) time.
   */
  if((r = tmap_fibheap_comparedata(h, key, data, x)) > 0) {
      /* XXX - bad code! */
      abort();
      tmap_fibheap_deleteel(h, x);

      x->tmap_fibheap_element_data = data;
      x->tmap_fibheap_element_key = key;

      tmap_fibheap_insertel(h, x);

      return odata;
  }

  x->tmap_fibheap_element_data = data;
  x->tmap_fibheap_element_key = key;

  /* because they are equal, we don't have to do anything */
  if(r == 0)
    return odata;

  y = x->tmap_fibheap_element_p;

  if(1 == h->tmap_fibheap_keys && okey == key)
    return odata;

  if(y != NULL && tmap_fibheap_compare(h, x, y) <= 0) {
      tmap_fibheap_cut(h, x, y);
      tmap_fibheap_cascading_cut(h, y);
  }

  /*
   * the = is so that the call from tmap_fibheap_delete will delete the proper
   * element.
   */
  if(tmap_fibheap_compare(h, x, h->tmap_fibheap_min) <= 0)
    h->tmap_fibheap_min = x;

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
tmap_fibheap_element_t *
tmap_fibheap_insert(tmap_fibheap_t *h, void *data)
{
  tmap_fibheap_element_t *x;

  if(1 == h->tmap_fibheap_keys) tmap_error("void heap required", Exit, OutOfRange);
  if((x = tmap_fibheap_element_newelem()) == NULL)
    return NULL;

  /* just insert on root list, and make sure it's not the new min */
  x->tmap_fibheap_element_data = data;

  tmap_fibheap_insertel(h, x);

  return x;
}

void *
tmap_fibheap_min(tmap_fibheap_t *h)
{
  if(h->tmap_fibheap_min == NULL)
    return NULL;
  return h->tmap_fibheap_min->tmap_fibheap_element_data;
}

void *
tmap_fibheap_extractmin(tmap_fibheap_t *h)
{
  tmap_fibheap_element_t *z;
  void *ret = NULL;

  if(h->tmap_fibheap_min != NULL) {
      z = tmap_fibheap_extractminel(h);
      ret = z->tmap_fibheap_element_data;
      tmap_fibheap_element_destroy(z);
  }

  return ret;
}

void *
tmap_fibheap_replacedata(tmap_fibheap_t *h, tmap_fibheap_element_t *x, void *data)
{
  return tmap_fibheap_replacekeydata(h, x, x->tmap_fibheap_element_key, data);
}

void *
tmap_fibheap_delete(tmap_fibheap_t *h, tmap_fibheap_element_t *x)
{
  void *k;

  k = x->tmap_fibheap_element_data;
  if(0 == h->tmap_fibheap_keys)
    tmap_fibheap_replacedata(h, x, h->tmap_fibheap_neginf);
  else
    tmap_fibheap_replacekey(h, x, INT_MIN);
  tmap_fibheap_extractmin(h);

  return k;
}

