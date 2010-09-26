/* The MIT License

   Copyright (c) 2008, by Attractive Chaos <attractivechaos@aol.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   */

#ifndef FMAP_SORT_H_
#define FMAP_SORT_H_

#include <stdlib.h>
#include <string.h>
#include "fmap_alloc.h"

/*! 
  Generic Sort Library
  */

/*! 
  stack structure for introsort
  @param  left   the left boundary
  @param  right   the right boundary
  @param  depth  the stack depth
  */
typedef struct {
    void *left;
    void *right;
    int depth;
} fmap_sort_isort_stack_t;

/*! 
  @param  type_t  the type of values [type]
  @param  a       the first variable
  @param  b       the second variable
*/
#define FMAP_SORT_SWAP(type_t, a, b) { register type_t t=(a); (a)=(b); (b)=t; }
/*! 
  initializes sort functions with the given name, type, and comparison function 
  @param  name       the name of sort functions [symbol]
  @param  type_t     the type of values [type]
  @param  __sort_lt  the comparison function
  */
#define FMAP_SORT_INIT(name, type_t, __sort_lt) \
  void fmap_sort_mergesort_##name(size_t n, type_t array[], type_t temp[]) \
{ \
  type_t *a2[2], *a, *b; \
  int curr, shift; \
  \
  a2[0] = array; \
  a2[1] = temp? temp : (type_t*)fmap_malloc(sizeof(type_t) * n, "a2[1]"); \
  for (curr = 0, shift = 0; (1ul<<shift) < n; ++shift) { \
      a = a2[curr]; b = a2[1-curr]; \
      if (shift == 0) { \
          type_t *p = b, *i, *eb = a + n; \
          for (i = a; i < eb; i += 2) { \
              if (i == eb - 1) *p++ = *i; \
              else { \
                  if (__sort_lt(*(i+1), *i)) { \
                      *p++ = *(i+1); *p++ = *i; \
                  } else { \
                      *p++ = *i; *p++ = *(i+1); \
                  } \
              } \
          } \
      } else { \
          size_t i, step = 1ul<<shift; \
          for (i = 0; i < n; i += step<<1) { \
              type_t *p, *j, *k, *ea, *eb; \
              if (n < i + step) { \
                  ea = a + n; eb = a; \
              } else { \
                  ea = a + i + step; \
                  eb = a + (n < i + (step<<1)? n : i + (step<<1)); \
              } \
              j = a + i; k = a + i + step; p = b + i; \
              while (j < ea && k < eb) { \
                  if (__sort_lt(*k, *j)) *p++ = *k++; \
                  else *p++ = *j++; \
              } \
              while (j < ea) *p++ = *j++; \
              while (k < eb) *p++ = *k++; \
          } \
      } \
      curr = 1 - curr; \
  } \
  if (curr == 1) { \
      type_t *p = a2[0], *i = a2[1], *eb = array + n; \
      for (; p < eb; ++i) *p++ = *i; \
  } \
  if (temp == 0) free(a2[1]); \
} \
void fmap_sort_heapadjust_##name(size_t i, size_t n, type_t l[]) \
{ \
  size_t k = i; \
  type_t tmp = l[i]; \
  while ((k = (k << 1) + 1) < n) { \
      if (k != n - 1 && __sort_lt(l[k], l[k+1])) ++k; \
      if (__sort_lt(l[k], tmp)) break; \
      l[i] = l[k]; i = k; \
  } \
  l[i] = tmp; \
} \
void fmap_sort_heapmake_##name(size_t lsize, type_t l[]) \
{ \
  size_t i; \
  for (i = (lsize >> 1) - 1; i != (size_t)(-1); --i) \
  fmap_sort_heapadjust_##name(i, lsize, l); \
} \
void fmap_sort_heapsort_##name(size_t lsize, type_t l[]) \
{ \
  size_t i; \
  for (i = lsize - 1; i > 0; --i) { \
      type_t tmp; \
      tmp = *l; *l = l[i]; l[i] = tmp; fmap_sort_heapadjust_##name(0, i, l); \
  } \
} \
inline void __fmap_sort_insertsort_##name(type_t *s, type_t *t) \
{ \
  type_t *i, *j, swap_tmp; \
  for (i = s + 1; i < t; ++i) \
  for (j = i; j > s && __sort_lt(*j, *(j-1)); --j) { \
      swap_tmp = *j; *j = *(j-1); *(j-1) = swap_tmp; \
  } \
} \
void fmap_sort_combsort_##name(size_t n, type_t a[]) \
{ \
  const double shrink_factor = 1.2473309501039786540366528676643; \
  int do_swap; \
  size_t gap = n; \
  type_t tmp, *i, *j; \
  do { \
      if (gap > 2) { \
          gap = (size_t)(gap / shrink_factor); \
          if (gap == 9 || gap == 10) gap = 11; \
      } \
      do_swap = 0; \
      for (i = a; i < a + n - gap; ++i) { \
          j = i + gap; \
          if (__sort_lt(*j, *i)) { \
              tmp = *i; *i = *j; *j = tmp; \
              do_swap = 1; \
          } \
      } \
  } while (do_swap || gap > 2); \
  if (gap != 1) __fmap_sort_insertsort_##name(a, a + n); \
} \
void fmap_sort_introsort_##name(size_t n, type_t a[]) \
{ \
  int d; \
  fmap_sort_isort_stack_t *top, *stack; \
  type_t rp, swap_tmp; \
  type_t *s, *t, *i, *j, *k; \
  \
  if (n < 1) return; \
  else if (n == 2) { \
      if (__sort_lt(a[1], a[0])) { swap_tmp = a[0]; a[0] = a[1]; a[1] = swap_tmp; } \
      return; \
  } \
  for (d = 2; 1ul<<d < n; ++d); \
  stack = (fmap_sort_isort_stack_t*)fmap_malloc(sizeof(fmap_sort_isort_stack_t) * ((sizeof(size_t)*d)+2), "stack"); \
  top = stack; s = a; t = a + (n-1); d <<= 1; \
  while (1) { \
      if (s < t) { \
          if (--d == 0) { \
              fmap_sort_combsort_##name(t - s + 1, s); \
              t = s; \
              continue; \
          } \
          i = s; j = t; k = i + ((j-i)>>1) + 1; \
          if (__sort_lt(*k, *i)) { \
              if (__sort_lt(*k, *j)) k = j; \
          } else k = __sort_lt(*j, *i)? i : j; \
          rp = *k; \
          if (k != t) { swap_tmp = *k; *k = *t; *t = swap_tmp; } \
          for (;;) { \
              do ++i; while (__sort_lt(*i, rp)); \
              do --j; while (i <= j && __sort_lt(rp, *j)); \
              if (j <= i) break; \
              swap_tmp = *i; *i = *j; *j = swap_tmp; \
          } \
          swap_tmp = *i; *i = *t; *t = swap_tmp; \
          if (i-s > t-i) { \
              if (i-s > 16) { top->left = s; top->right = i-1; top->depth = d; ++top; } \
              s = t-i > 16? i+1 : t; \
          } else { \
              if (t-i > 16) { top->left = i+1; top->right = t; top->depth = d; ++top; } \
              t = i-s > 16? i-1 : s; \
          } \
      } else { \
          if (top == stack) { \
              free(stack); \
              __fmap_sort_insertsort_##name(a, a+n); \
              return; \
          } else { --top; s = (type_t*)top->left; t = (type_t*)top->right; d = top->depth; } \
      } \
  } \
} \
/* This function is adapted from: http://ndevilla.free.fr/median/ */ \
/* 0 <= kk < n */ \
type_t fmap_sort_small_##name(size_t n, type_t arr[], size_t kk) \
{ \
  type_t *low, *high, *k, *ll, *hh, *mid; \
  low = arr; high = arr + n - 1; k = arr + kk; \
  for (;;) { \
      if (high <= low) return *k; \
      if (high == low + 1) { \
          if (__sort_lt(*high, *low)) FMAP_SORT_SWAP(type_t, *low, *high); \
          return *k; \
      } \
      mid = low + (high - low) / 2; \
      if (__sort_lt(*high, *mid)) FMAP_SORT_SWAP(type_t, *mid, *high); \
      if (__sort_lt(*high, *low)) FMAP_SORT_SWAP(type_t, *low, *high); \
      if (__sort_lt(*low, *mid)) FMAP_SORT_SWAP(type_t, *mid, *low); \
      FMAP_SORT_SWAP(type_t, *mid, *(low+1)); \
      ll = low + 1; hh = high; \
      for (;;) { \
          do ++ll; while (__sort_lt(*ll, *low)); \
          do --hh; while (__sort_lt(*low, *hh)); \
          if (hh < ll) break; \
          FMAP_SORT_SWAP(type_t, *ll, *hh); \
      } \
      FMAP_SORT_SWAP(type_t, *low, *hh); \
      if (hh <= k) low = ll; \
      if (hh >= k) high = hh - 1; \
  } \
}

/*! 
  performs mergesort on the given array
  @param  name  the name of the sort functions [symbol] 
  @param  n     the size of the array
  @param  a     the array of elements to be sorted
  @param  t     a temporary array of elements of length n, or NULL
  */
#define fmap_sort_mergesort(name, n, a, t) fmap_sort_mergesort_##name(n, a, t)
/*! 
  performs introsort on the given array
  @param  name  the name of the sort functions [symbol] 
  @param  n     the size of the array
  @param  a     the array of elements to be sorted
  */
#define fmap_sort_introsort(name, n, a) fmap_sort_introsort_##name(n, a)
/*! 
  performs combsort on the given array
  @param  name  the name of the sort functions [symbol] 
  @param  n     the size of the array
  @param  a     the array of elements to be sorted
  */
#define fmap_sort_combsort(name, n, a) fmap_sort_combsort_##name(n, a)
/*! 
  performs heapsort on the given array
  @param  name  the name of the sort functions [symbol] 
  @param  n     the size of the array
  @param  a     the array of elements to be sorted
  */
#define fmap_sort_heapsort(name, n, a) fmap_sort_heapsort_##name(n, a)
#define fmap_sort_heapmake(name, n, a) fmap_sort_heapmake_##name(n, a)
#define fmap_sort_heapadjust(name, i, n, a) fmap_sort_heapadjust_##name(i, n, a)
/*! 
  performs small sorton the given array
  @param  name  the name of the sort functions [symbol] 
  @param  n     the size of the array
  @param  a     the array of elements to be sorted
  @param  k     a midpoint to start this algorithm
  @details      adapted from http://ndevilla.free.fr/median/
  */
#define fmap_sort_small(name, n, a, k) fmap_sort_small_##name(n, a, k)

/*! 
  a generic value-based less-than comparison function
  @param  a  the first value to compare
  @param  b  the second value to compare
  @return    1 if true, 0 otherwise
  */
#define fmap_sort_lt_generic(a, b) ((a) < (b))
/*! 
  a generic string-based comparison function
  @param  a  the first string to compare
  @param  b  the second string to compare
  @return    1 if true, 0 otherwise
  @details   this uses strcmp
  */
#define fmap_sort_lt_str(a, b) (strcmp((a), (b)) < 0)

typedef const char *ksstr_t;

/*! 
  initializes sort functions with the given type
  @param  type_t  the type of values [type]
  @details        this will define functions named by the type and using the generic value-based comparison function.
  */
#define FMAP_SORT_INIT_GENERIC(type_t) FMAP_SORT_INIT(type_t, type_t, fmap_sort_lt_generic)

/*! 
  initializes string comparison sort functions 
  @details  this will define functions with the name "str" and using the generic string-based comparison function.
  */
#define FMAP_SORT_INIT_STR FMAP_SORT_INIT(str, ksstr_t, fmap_sort_lt_str)

#endif
