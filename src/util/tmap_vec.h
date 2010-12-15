/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
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

#ifndef TMAP_VEC_H
#define TMAP_VEC_H

#include <stdlib.h>
#include "tmap_alloc.h"

/*! 
  A Vector Library
  */

#ifdef tmap_roundup32
#define tmap_vec_roundup32 tmap_roundup32
#else
/*! 
  rounds up to the nearest power of two integer
  @param  x  the integer to round up
  @return    the smallest integer greater than x that is a power of two 
  */
#define tmap_vec_roundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/*!   tmap_vec_t(type)
  @param  type  the type of values [type]
  @details      defines an anonymous struct representing a vector of values of the given type  
*/
#define tmap_vec_t(type) struct { size_t n, m; type *a; }
/*! 
  @param  v  the vector to initialize
*/
#define tmap_vec_init(v) ((v).n = (v).m = 0, (v).a = 0)
/*! 
  @param  v  the vector to destroy 
*/
#define tmap_vec_destroy(v) free((v).a)
/*! 
  sets the element at index i
  @param  v  the element to set
  @param  i  the zero-based index
  @return    the set element
*/
#define tmap_vec_A(v, i) ((v).a[(i)])
/*! 
  pop an element from the vector
  @param  v  the vector 
  @return    the popped element 
*/
#define tmap_vec_pop(v) ((v).a[--(v).n])
/*! 
  get the vector size
  @param  v  the vector 
  @return    the size of the vector
*/
#define tmap_vec_size(v) ((v).n)
/*! 
  get the vectors memory capacity
  @param  v  the vector 
  @return    the memory capacity
*/
#define tmap_vec_max(v) ((v).m)
/*! 
  resizes the given vector
  @param  type  the type of values [type]
  @param  v     the vector 
  @param  s     the new memory size
*/
#define tmap_vec_resize(type, v, s)  ((v).m = (s), (v).a = (type*)tmap_realloc((v).a, sizeof(type) * (v).m, "(v).a"))
/*! 
  copies vector v0 into vector v1
  @param  type  the type of values [type]
  @param  v1    the destination vector
  @param  v0    the source vector
*/
#define tmap_vec_copy(type, v1, v0) do { \
  if ((v1).m < (v0).n) tmap_vec_resize(type, v1, (v0).n); \
  (v1).n = (v0).n; \
  memcpy((v1).a, (v0).a, sizeof(type) * (v0).n); \
} while (0) \
/*! 
  adds the give element to the vector
  @param  type  the type of values [type]
  @param  v     the vector
  @param  x     the element to add
  */
#define tmap_vec_push(type, v, x) do { \
    if ((v).n == (v).m) { \
        (v).m = (v).m? (v).m<<1 : 2; \
        (v).a = (type*)tmap_realloc((v).a, sizeof(type) * (v).m, "(v).a"); \
    } \
    (v).a[(v).n++] = (x); \
} while (0)

// Undocumented
#define tmap_vec_pushp(type, v) (((v).n == (v).m)? \
                                 ((v).m = ((v).m? (v).m<<1 : 2), \
                                  (v).a = (type*)realloc((v).a, sizeof(type) * (v).m, "(v).a"), 0) \
                                 : 0), ((v).a + ((v).n++))

// Undocumented
#define tmap_vec_a(type, v, i) ((v).m <= (size_t)(i)? \
                                ((v).m = (v).n = (i) + 1, tmap_vec_roundup32((v).m), \
                                 (v).a = (type*)tmap_realloc((v).a, sizeof(type) * (v).m, "(v).a"), 0) \
                                : (v).n <= (size_t)(i)? (v).n = (i) \
                                : 0), (v).a[(i)]

#endif
