/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/* The MIT License

   Copyright (c) 2008, 2009 by attractor <attractor@live.co.uk>

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


#ifndef TMAP_HASH_H
#define TMAP_HASH_H

/*! 
  Generic Hash Libary
  */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tmap_alloc.h"

/* compipler specific configuration */

typedef uint32_t tmap_hash_int32_t;
typedef uint64_t tmap_hash_int64_t;
typedef tmap_hash_int32_t tmap_hash_int_t;
typedef tmap_hash_int_t tmap_hash_iter_t;

#define TMAP_HASH_PRIME_SIZE 32
static const tmap_hash_int32_t __tmap_hash_prime_list[TMAP_HASH_PRIME_SIZE] =
{
  0ul,          3ul,          11ul,         23ul,         53ul,
  97ul,         193ul,        389ul,        769ul,        1543ul,
  3079ul,       6151ul,       12289ul,      24593ul,      49157ul,
  98317ul,      196613ul,     393241ul,     786433ul,     1572869ul,
  3145739ul,    6291469ul,    12582917ul,   25165843ul,   50331653ul,
  100663319ul,  201326611ul,  402653189ul,  805306457ul,  1610612741ul,
  3221225473ul, 4294967291ul
};

#define __tmap_hash_isempty(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&2)
#define __tmap_hash_isdel(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&1)
#define __tmap_hash_iseither(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&3)
#define __tmap_hash_set_isdel_false(flag, i) (flag[i>>4]&=~(1ul<<((i&0xfU)<<1)))
#define __tmap_hash_set_isempty_false(flag, i) (flag[i>>4]&=~(2ul<<((i&0xfU)<<1)))
#define __tmap_hash_set_isboth_false(flag, i) (flag[i>>4]&=~(3ul<<((i&0xfU)<<1)))
#define __tmap_hash_set_isdel_true(flag, i) (flag[i>>4]|=1ul<<((i&0xfU)<<1))

static const double __tmap_hash_HASH_UPPER = 0.77;

#define TMAP_HASH_INIT(name, tmap_hashkey_t, tmap_hashval_t, tmap_hash_is_map, __hash_func, __hash_equal) \
  typedef struct { \
      tmap_hash_int_t n_buckets, size, n_occupied, upper_bound; \
      tmap_hash_int32_t *flags; \
      tmap_hashkey_t *keys; \
      tmap_hashval_t *vals; \
  } tmap_hash_##name##_t; \
static inline tmap_hash_##name##_t *tmap_hash_init_##name() { \
    return (tmap_hash_##name##_t*)tmap_calloc(1, sizeof(tmap_hash_##name##_t), "tmap_hash"); \
} \
static inline void tmap_hash_destroy_##name(tmap_hash_##name##_t *h) \
{ \
  if (h) { \
      free(h->keys); free(h->flags); \
      free(h->vals); \
      free(h); \
  } \
} \
static inline void tmap_hash_clear_##name(tmap_hash_##name##_t *h) \
{ \
  if (h && NULL != h->flags) { \
      memset(h->flags, 0xaa, ((h->n_buckets>>4) + 1) * sizeof(tmap_hash_int32_t)); \
      h->size = h->n_occupied = 0; \
  } \
} \
static inline tmap_hash_int_t tmap_hash_get_##name(const tmap_hash_##name##_t *h, tmap_hashkey_t key) \
{ \
  if (h->n_buckets) { \
      tmap_hash_int_t inc, k, i, last; \
      k = __hash_func(key); i = k % h->n_buckets; \
      inc = 1 + k % (h->n_buckets - 1); last = i; \
      while (!__tmap_hash_isempty(h->flags, i) && (__tmap_hash_isdel(h->flags, i) || !__hash_equal(h->keys[i], key))) { \
          if (i + inc >= h->n_buckets) i = i + inc - h->n_buckets; \
          else i += inc; \
          if (i == last) return h->n_buckets; \
      } \
      return __tmap_hash_iseither(h->flags, i)? h->n_buckets : i; \
  } else return 0; \
} \
\
static inline void \
tmap_hash_resize_##name(tmap_hash_##name##_t *h, tmap_hash_int_t new_n_buckets) \
{ \
  tmap_hash_int32_t *new_flags = 0; \
  tmap_hash_int_t j = 1; \
    { \
      tmap_hash_int_t t = TMAP_HASH_PRIME_SIZE - 1; \
      while (__tmap_hash_prime_list[t] > new_n_buckets) --t; \
      new_n_buckets = __tmap_hash_prime_list[t+1]; \
      if (h->size >= (tmap_hash_int_t)(new_n_buckets * __tmap_hash_HASH_UPPER + 0.5)) j = 0; \
      else { \
          new_flags = (tmap_hash_int32_t*)tmap_malloc(((new_n_buckets>>4) + 1) * sizeof(tmap_hash_int32_t), "new_flags"); \
          memset(new_flags, 0xaa, ((new_n_buckets>>4) + 1) * sizeof(tmap_hash_int32_t)); \
          if (h->n_buckets < new_n_buckets) { \
              h->keys = (tmap_hashkey_t*)tmap_realloc(h->keys, new_n_buckets * sizeof(tmap_hashkey_t), "h->keys"); \
              if (tmap_hash_is_map) \
              h->vals = (tmap_hashval_t*)tmap_realloc(h->vals, new_n_buckets * sizeof(tmap_hashval_t), "h->vals"); \
          } \
      } \
    } \
  if (j) { \
      for (j = 0; j != h->n_buckets; ++j) { \
          if (__tmap_hash_iseither(h->flags, j) == 0) { \
              tmap_hashkey_t key = h->keys[j]; \
              tmap_hashval_t val; \
              if (tmap_hash_is_map) val = h->vals[j]; \
              __tmap_hash_set_isdel_true(h->flags, j); \
              while (1) { \
                  tmap_hash_int_t inc, k, i; \
                  k = __hash_func(key); \
                  i = k % new_n_buckets; \
                  inc = 1 + k % (new_n_buckets - 1); \
                  while (!__tmap_hash_isempty(new_flags, i)) { \
                      if (i + inc >= new_n_buckets) i = i + inc - new_n_buckets; \
                      else i += inc; \
                  } \
                  __tmap_hash_set_isempty_false(new_flags, i); \
                  if (i < h->n_buckets && __tmap_hash_iseither(h->flags, i) == 0) { \
                        { tmap_hashkey_t tmp = h->keys[i]; h->keys[i] = key; key = tmp; } \
                      if (tmap_hash_is_map) { tmap_hashval_t tmp = h->vals[i]; h->vals[i] = val; val = tmp; } \
                      __tmap_hash_set_isdel_true(h->flags, i); \
                  } else { \
                      h->keys[i] = key; \
                      if (tmap_hash_is_map) h->vals[i] = val; \
                      break; \
                  } \
              } \
          } \
      } \
      if (h->n_buckets > new_n_buckets) { \
          h->keys = (tmap_hashkey_t*)tmap_realloc(h->keys, new_n_buckets * sizeof(tmap_hashkey_t), "h->keys"); \
          if (tmap_hash_is_map) \
          h->vals = (tmap_hashval_t*)tmap_realloc(h->vals, new_n_buckets * sizeof(tmap_hashval_t), "h->vals"); \
      } \
      free(h->flags); \
      h->flags = new_flags; \
      h->n_buckets = new_n_buckets; \
      h->n_occupied = h->size; \
      h->upper_bound = (tmap_hash_int_t)(h->n_buckets * __tmap_hash_HASH_UPPER + 0.5); \
  } \
} \
static inline tmap_hash_int_t tmap_hash_put_##name(tmap_hash_##name##_t *h, tmap_hashkey_t key, int *ret) \
{ \
  tmap_hash_int_t x; \
  if (h->n_occupied >= h->upper_bound) { \
      if (h->n_buckets > (h->size<<1)) tmap_hash_resize_##name(h, h->n_buckets - 1); \
      else tmap_hash_resize_##name(h, h->n_buckets + 1); \
  } \
    { \
      tmap_hash_int_t inc, k, i, site, last; \
      x = site = h->n_buckets; k = __hash_func(key); i = k % h->n_buckets; \
      if (__tmap_hash_isempty(h->flags, i)) x = i; \
      else { \
          inc = 1 + k % (h->n_buckets - 1); last = i; \
          while (!__tmap_hash_isempty(h->flags, i) && (__tmap_hash_isdel(h->flags, i) || !__hash_equal(h->keys[i], key))) { \
              if (__tmap_hash_isdel(h->flags, i)) site = i; \
              if (i + inc >= h->n_buckets) i = i + inc - h->n_buckets; \
              else i += inc; \
              if (i == last) { x = site; break; } \
          } \
          if (x == h->n_buckets) { \
              if (__tmap_hash_isempty(h->flags, i) && site != h->n_buckets) x = site; \
              else x = i; \
          } \
      } \
    } \
  if (__tmap_hash_isempty(h->flags, x)) { \
      h->keys[x] = key; \
      __tmap_hash_set_isboth_false(h->flags, x); \
      ++h->size; ++h->n_occupied; \
      *ret = 1; \
  } else if (__tmap_hash_isdel(h->flags, x)) { \
      h->keys[x] = key; \
      __tmap_hash_set_isboth_false(h->flags, x); \
      ++h->size; \
      *ret = 2; \
  } else *ret = 0; \
  return x; \
} \
static inline void tmap_hash_del_##name(tmap_hash_##name##_t *h, tmap_hash_int_t x) \
{ \
  if (x != h->n_buckets && !__tmap_hash_iseither(h->flags, x)) { \
      __tmap_hash_set_isdel_true(h->flags, x); \
      --h->size; \
  } \
}

/* --- BEGIN OF HASH FUNCTIONS --- */

/*! 
  Integer hash function
  @param  key   The integer [tmap_hash_int32_t]
  @return       The hash value [tmap_hash_int_t]
  */
#define tmap_hash_int_hash_func(key) (tmap_hash_int32_t)(key)
/*! 
  Integer comparison function
  */
#define tmap_hash_int_hash_equal(a, b) ((a) == (b))
/*! 
  64-bit integer hash function
  @param  key   The integer [tmap_hash_int64_t]
  @return       The hash value [tmap_hash_int_t]
  */
#define tmap_hash_int64_hash_func(key) (tmap_hash_int32_t)((key)>>33^(key)^(key)<<11)
/*! 
  64-bit integer comparison function
  */
#define tmap_hash_int64_hash_equal(a, b) ((a) == (b))
// Undocumented
static inline tmap_hash_int_t __tmap_hash_X31_hash_string(const char *s)
{
  tmap_hash_int_t h = *s;
  if (h) for (++s ; *s; ++s) h = (h << 5) - h + *s;
  return h;
}
/*! 
  Another interface to const char* hash function
  @param  key   Pointer to a null terminated string [const char*]
  @return       The hash value [tmap_hash_int_t]
  */
#define tmap_hash_str_hash_func(key) __tmap_hash_X31_hash_string(key)
/*! 
  Const char* comparison function
  */
#define tmap_hash_str_hash_equal(a, b) (strcmp(a, b) == 0)

/* --- END OF HASH FUNCTIONS --- */

/* Other necessary defines... */

/*!
  @param  name  Name of the hash table [symbol]
  */
#define tmap_hash_t(name) tmap_hash_##name##_t

/*! 
  Initiate a hash table.
  @param  name  Name of the hash table [symbol]
  @return       Pointer to the hash table [tmap_hash_t(name)*]
  */
#define tmap_hash_init(name) tmap_hash_init_##name()

/*! 
  Destroy a hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  */
#define tmap_hash_destroy(name, h) tmap_hash_destroy_##name(h)

/*! 
  Reset a hash table without deallocating memory.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  */
#define tmap_hash_clear(name, h) tmap_hash_clear_##name(h)

/*! 
  Resize a hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @param  s     New size [tmap_hash_int_t]
  */
#define tmap_hash_resize(name, h, s) tmap_hash_resize_##name(h, s)

/*! 
  Insert a key to the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @param  k     Key [type of keys]
  @param  r     Extra return code: 0 if the key is present in the hash table;
  1 if the bucket is empty (never used); 2 if the element in
  the bucket has been deleted [int*]
  @return       Iterator to the inserted element [tmap_hash_int_t]
  */
#define tmap_hash_put(name, h, k, r) tmap_hash_put_##name(h, k, r)

/*! 
  Retrieve a key from the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @param  k     Key [type of keys]
  @return       Iterator to the found element, or tmap_hash_end(h) is the element is absent [tmap_hash_int_t]
  */
#define tmap_hash_get(name, h, k) tmap_hash_get_##name(h, k)

/*! 
  Remove a key from the hash table.
  @param  name  Name of the hash table [symbol]
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @param  k     Iterator to the element to be deleted [tmap_hash_int_t]
  */
#define tmap_hash_del(name, h, k) tmap_hash_del_##name(h, k)


/*! 
  Test whether a bucket contains data.
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @param  x     Iterator to the bucket [tmap_hash_int_t]
  @return       1 if containing data; 0 otherwise [int]
  */
#define tmap_hash_exist(h, x) (!__tmap_hash_iseither((h)->flags, (x)))

/*! 
  Get key given an iterator
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @param  x     Iterator to the bucket [tmap_hash_int_t]
  @return       Key [type of keys]
  */
#define tmap_hash_key(h, x) ((h)->keys[x])

/*! 
  Get value given an iterator
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @param  x     Iterator to the bucket [tmap_hash_int_t]
  @return       Value [type of values]
  @details      For hash sets, calling this results in segfault.
  */
#define tmap_hash_val(h, x) ((h)->vals[x])

/*! 
  Alias of tmap_hash_val()
  */
#define tmap_hash_value(h, x) ((h)->vals[x])

/*! 
  Get the start iterator
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @return       The start iterator [tmap_hash_int_t]
  */
#define tmap_hash_begin(h) (tmap_hash_int_t)(0)

/*! 
  Get the end iterator
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @return       The end iterator [tmap_hash_int_t]
  */
#define tmap_hash_end(h) ((h)->n_buckets)

/*! 
  Get the number of elements in the hash table
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @return       Number of elements in the hash table [tmap_hash_int_t]
  */
#define tmap_hash_size(h) ((h)->size)

/*! 
  Get the number of buckets in the hash table
  @param  h     Pointer to the hash table [tmap_hash_t(name)*]
  @return       Number of buckets in the hash table [tmap_hash_int_t]
  */
#define tmap_hash_n_buckets(h) ((h)->n_buckets)

/* More conenient interfaces */

/*! 
  Instantiate a hash set containing integer keys
  @param  name  Name of the hash table [symbol]
  */
#define TMAP_HASH_SET_INIT_INT(name) \
  TMAP_HASH_INIT(name, tmap_hash_int32_t, char, 0, tmap_hash_int_hash_func, tmap_hash_int_hash_equal)

/*! 
  Instantiate a hash map containing integer keys
  @param  name  Name of the hash table [symbol]
  @param  tmap_hashval_t  Type of values [type]
  */
#define TMAP_HASH_MAP_INIT_INT(name, tmap_hashval_t) \
  TMAP_HASH_INIT(name, tmap_hash_int32_t, tmap_hashval_t, 1, tmap_hash_int_hash_func, tmap_hash_int_hash_equal)

/*! 
  Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
  */
#define TMAP_HASH_SET_INIT_INT64(name) \
  TMAP_HASH_INIT(name, tmap_hash_int64_t, char, 0, tmap_hash_int64_hash_func, tmap_hash_int64_hash_equal)

/*! 
  Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
  @param  tmap_hashval_t  Type of values [type]
  */
#define TMAP_HASH_MAP_INIT_INT64(name, tmap_hashval_t) \
  TMAP_HASH_INIT(name, tmap_hash_int64_t, tmap_hashval_t, 1, tmap_hash_int64_hash_func, tmap_hash_int64_hash_equal)

typedef const char *tmap_hash_cstr_t;
/*! 
  Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
  */
#define TMAP_HASH_SET_INIT_STR(name) \
  TMAP_HASH_INIT(name, tmap_hash_cstr_t, char, 0, tmap_hash_str_hash_func, tmap_hash_str_hash_equal)

/*! 
  Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
  @param  tmap_hashval_t  Type of values [type]
  */
#define TMAP_HASH_MAP_INIT_STR(name, tmap_hashval_t) \
  TMAP_HASH_INIT(name, tmap_hash_cstr_t, tmap_hashval_t, 1, tmap_hash_str_hash_func, tmap_hash_str_hash_equal)

#endif 
