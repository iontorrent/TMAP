/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/* The MIT License

   Copyright (c) 2008, 2009, 2011 by Attractive Chaos <attractor@live.co.uk>

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

/*
  An example:

#include "tmap_hash.h"
TMAP_HASH_MAP_INIT_INT(32, char)
int main() {
	int ret, is_missing;
	tmap_iter_t k;
	tmap_hash_t(32) *h = tmap_hash_init(32);
	k = tmap_hash_put(32, h, 5, &ret);
	if (!ret) tmap_hash_del(32, h, k);
	tmap_hash_value(h, k) = 10;
	k = tmap_hash_get(32, h, 10);
	is_missing = (k == tmap_hash_end(h));
	k = tmap_hash_get(32, h, 5);
	tmap_hash_del(32, h, k);
	for (k = tmap_hash_begin(h); k != tmap_hash_end(h); ++k)
		if (tmap_hash_exist(h, k)) tmap_hash_value(h, k) = 1;
	tmap_hash_destroy(32, h);
	return 0;
}
*/

/*
  2011-09-16 (0.2.6):

	* The capacity is a power of 2. This seems to dramatically improve the
	  speed for simple keys. Thank Zilong Tan for the suggestion. Reference:

	   - http://code.google.com/p/ulib/
	   - http://nothings.org/computer/judy/

	* Allow to optionally use linear probing which usually has better
	  performance for random input. Double hashing is still the default as it
	  is more robust to certain non-random input.

	* Added Wang's integer hash function (not used by default). This hash
	  function is more robust to certain non-random input.

  2011-02-14 (0.2.5):

    * Allow to declare global functions.

  2009-09-26 (0.2.4):

    * Improve portability

  2008-09-19 (0.2.3):

	* Corrected the example
	* Improved interfaces

  2008-09-11 (0.2.2):

	* Improved speed a little in tmap_hash_put()

  2008-09-10 (0.2.1):

	* Added tmap_hash_clear()
	* Fixed a compiling error

  2008-09-02 (0.2.0):

	* Changed to token concatenation which increases flexibility.

  2008-08-31 (0.1.2):

	* Fixed a bug in tmap_hash_get(), which has not been tested previously.

  2008-08-31 (0.1.1):

	* Added destructor
*/

#ifndef TMAP_HASH_H
#define TMAP_HASH_H

/*! 
  Generic Hash Libary
  */

#define TMAP_HASH_VERSION_TMAP_HASH_H "0.2.6"

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tmap_alloc.h"
#include "tmap_definitions.h"

/* compiler specific configuration */

typedef uint32_t tmap_hash_int32_t;
typedef uint64_t tmap_hash_int64_t;
typedef tmap_hash_int32_t tmap_hash_int_t;
typedef tmap_hash_int_t tmap_hash_iter_t;

#define __tmap_hash_isempty(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&2)
#define __tmap_hash_isdel(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&1)
#define __tmap_hash_iseither(flag, i) ((flag[i>>4]>>((i&0xfU)<<1))&3)
#define __tmap_hash_set_isdel_false(flag, i) (flag[i>>4]&=~(1ul<<((i&0xfU)<<1)))
#define __tmap_hash_set_isempty_false(flag, i) (flag[i>>4]&=~(2ul<<((i&0xfU)<<1)))
#define __tmap_hash_set_isboth_false(flag, i) (flag[i>>4]&=~(3ul<<((i&0xfU)<<1)))
#define __tmap_hash_set_isdel_true(flag, i) (flag[i>>4]|=1ul<<((i&0xfU)<<1))

#ifdef TMAP_HASH_LINEAR
#define __tmap_hash_inc(k, m) 1
#else
#define __tmap_hash_inc(k, m) (((k)>>3 ^ (k)<<3) | 1) & (m)
#endif

#define __tmap_hash_fsize(m) ((m) < 16? 1 : (m)>>4)

static const double __tmap_hash_HASH_UPPER = 0.77;

#define TMAP_HASH_DECLARE(name, tmap_key_t, tmap_val_t)		 					\
	typedef struct {													\
		tmap_hash_int_t n_buckets, size, n_occupied, upper_bound;				\
		int32_t *flags;												\
		tmap_key_t *keys;													\
		tmap_val_t *vals;													\
	} tmap_hash_##name##_t;													\
	extern tmap_hash_##name##_t *tmap_hash_init_##name();								\
	extern void tmap_hash_destroy_##name(tmap_hash_##name##_t *h);					\
	extern void tmap_hash_clear_##name(tmap_hash_##name##_t *h);						\
	extern tmap_hash_int_t tmap_hash_get_##name(const tmap_hash_##name##_t *h, tmap_key_t key); 	\
	extern void tmap_hash_resize_##name(tmap_hash_##name##_t *h, tmap_hash_int_t new_n_buckets); \
	extern tmap_hash_int_t tmap_hash_put_##name(tmap_hash_##name##_t *h, tmap_key_t key, int *ret); \
	extern void tmap_hash_del_##name(tmap_hash_##name##_t *h, tmap_hash_int_t x);

#define TMAP_HASH_INIT2(name, SCOPE, tmap_key_t, tmap_val_t, tmap_hash_is_map, __hash_func, __hash_equal) \
	typedef struct {													\
		tmap_hash_int_t n_buckets, size, n_occupied, upper_bound;				\
		int32_t *flags;												\
		tmap_key_t *keys;													\
		tmap_val_t *vals;													\
	} tmap_hash_##name##_t;													\
	SCOPE tmap_hash_##name##_t *tmap_hash_init_##name() {								\
		return (tmap_hash_##name##_t*)tmap_calloc(1, sizeof(tmap_hash_##name##_t), "return");		\
	}																	\
	SCOPE void tmap_hash_destroy_##name(tmap_hash_##name##_t *h)						\
	{																	\
		if (h) {														\
			free(h->keys); free(h->flags);								\
			free(h->vals);												\
			free(h);													\
		}																\
	}																	\
	SCOPE void tmap_hash_clear_##name(tmap_hash_##name##_t *h)						\
	{																	\
		if (h && h->flags) {											\
			memset(h->flags, 0xaa, __tmap_hash_fsize(h->n_buckets) * sizeof(int32_t)); \
			h->size = h->n_occupied = 0;								\
		}																\
	}																	\
	SCOPE tmap_hash_int_t tmap_hash_get_##name(const tmap_hash_##name##_t *h, tmap_key_t key) 	\
	{																	\
		if (h->n_buckets) {												\
			tmap_hash_int_t inc, k, i, last, mask;								\
			mask = h->n_buckets - 1;									\
			k = __hash_func(key); i = k & mask;							\
			inc = __tmap_hash_inc(k, mask); last = i; /* inc==1 for linear probing */ \
			while (!__tmap_hash_isempty(h->flags, i) && (__tmap_hash_isdel(h->flags, i) || !__hash_equal(h->keys[i], key))) { \
				i = (i + inc) & mask; 									\
				if (i == last) return h->n_buckets;						\
			}															\
			return __tmap_hash_iseither(h->flags, i)? h->n_buckets : i;		\
		} else return 0;												\
	}																	\
	SCOPE void tmap_hash_resize_##name(tmap_hash_##name##_t *h, tmap_hash_int_t new_n_buckets) \
	{ /* This function uses 0.25*n_bucktes bytes of working space instead of [sizeof(key_t+val_t)+.25]*n_buckets. */ \
		int32_t *new_flags = 0;										\
		tmap_hash_int_t j = 1;													\
		{																\
			tmap_roundup32(new_n_buckets); 									\
			if (new_n_buckets < 4) new_n_buckets = 4;					\
			if (h->size >= (tmap_hash_int_t)(new_n_buckets * __tmap_hash_HASH_UPPER + 0.5)) j = 0;	/* requested size is too small */ \
			else { /* hash table size to be changed (shrink or expand); rehash */ \
				new_flags = (int32_t*)tmap_malloc(__tmap_hash_fsize(new_n_buckets) * sizeof(int32_t), "new_flags");	\
				memset(new_flags, 0xaa, __tmap_hash_fsize(new_n_buckets) * sizeof(int32_t)); \
				if (h->n_buckets < new_n_buckets) {	/* expand */		\
					h->keys = (tmap_key_t*)tmap_realloc(h->keys, new_n_buckets * sizeof(tmap_key_t), "h->keys"); \
					if (tmap_hash_is_map) h->vals = (tmap_val_t*)tmap_realloc(h->vals, new_n_buckets * sizeof(tmap_val_t), "h->vals"); \
				} /* otherwise shrink */								\
			}															\
		}																\
		if (j) { /* rehashing is needed */								\
			for (j = 0; j != h->n_buckets; ++j) {						\
				if (__tmap_hash_iseither(h->flags, j) == 0) {					\
					tmap_key_t key = h->keys[j];							\
					tmap_val_t val;										\
					tmap_hash_int_t new_mask;									\
					new_mask = new_n_buckets - 1; 						\
					if (tmap_hash_is_map) val = h->vals[j];					\
					__tmap_hash_set_isdel_true(h->flags, j);					\
					while (1) { /* kick-out process; sort of like in Cuckoo hashing */ \
						tmap_hash_int_t inc, k, i;								\
						k = __hash_func(key);							\
						i = k & new_mask;								\
						inc = __tmap_hash_inc(k, new_mask);					\
						while (!__tmap_hash_isempty(new_flags, i)) i = (i + inc) & new_mask; \
						__tmap_hash_set_isempty_false(new_flags, i);			\
						if (i < h->n_buckets && __tmap_hash_iseither(h->flags, i) == 0) { /* kick out the existing element */ \
							{ tmap_key_t tmp = h->keys[i]; h->keys[i] = key; key = tmp; } \
							if (tmap_hash_is_map) { tmap_val_t tmp = h->vals[i]; h->vals[i] = val; val = tmp; } \
							__tmap_hash_set_isdel_true(h->flags, i); /* mark it as deleted in the old hash table */ \
						} else { /* write the element and jump out of the loop */ \
							h->keys[i] = key;							\
							if (tmap_hash_is_map) h->vals[i] = val;			\
							break;										\
						}												\
					}													\
				}														\
			}															\
			if (h->n_buckets > new_n_buckets) { /* shrink the hash table */ \
				h->keys = (tmap_key_t*)tmap_realloc(h->keys, new_n_buckets * sizeof(tmap_key_t), "h->keys"); \
				if (tmap_hash_is_map) h->vals = (tmap_val_t*)tmap_realloc(h->vals, new_n_buckets * sizeof(tmap_val_t), "h->vals"); \
			}															\
			free(h->flags); /* free the working space */				\
			h->flags = new_flags;										\
			h->n_buckets = new_n_buckets;								\
			h->n_occupied = h->size;									\
			h->upper_bound = (tmap_hash_int_t)(h->n_buckets * __tmap_hash_HASH_UPPER + 0.5); \
		}																\
	}																	\
	SCOPE tmap_hash_int_t tmap_hash_put_##name(tmap_hash_##name##_t *h, tmap_key_t key, int *ret) \
	{																	\
		tmap_hash_int_t x;														\
		if (h->n_occupied >= h->upper_bound) { /* update the hash table */ \
			if (h->n_buckets > (h->size<<1)) tmap_hash_resize_##name(h, h->n_buckets - 1); /* clear "deleted" elements */ \
			else tmap_hash_resize_##name(h, h->n_buckets + 1); /* expand the hash table */ \
		} /* TODO: to implement automatically shrinking; resize() already support shrinking */ \
		{																\
			tmap_hash_int_t inc, k, i, site, last, mask = h->n_buckets - 1;		\
			x = site = h->n_buckets; k = __hash_func(key); i = k & mask; \
			if (__tmap_hash_isempty(h->flags, i)) x = i; /* for speed up */	\
			else {														\
				inc = __tmap_hash_inc(k, mask); last = i;						\
				while (!__tmap_hash_isempty(h->flags, i) && (__tmap_hash_isdel(h->flags, i) || !__hash_equal(h->keys[i], key))) { \
					if (__tmap_hash_isdel(h->flags, i)) site = i;				\
					i = (i + inc) & mask; 								\
					if (i == last) { x = site; break; }					\
				}														\
				if (x == h->n_buckets) {								\
					if (__tmap_hash_isempty(h->flags, i) && site != h->n_buckets) x = site; \
					else x = i;											\
				}														\
			}															\
		}																\
		if (__tmap_hash_isempty(h->flags, x)) { /* not present at all */		\
			h->keys[x] = key;											\
			__tmap_hash_set_isboth_false(h->flags, x);							\
			++h->size; ++h->n_occupied;									\
			*ret = 1;													\
		} else if (__tmap_hash_isdel(h->flags, x)) { /* deleted */				\
			h->keys[x] = key;											\
			__tmap_hash_set_isboth_false(h->flags, x);							\
			++h->size;													\
			*ret = 2;													\
		} else *ret = 0; /* Don't touch h->keys[x] if present and not deleted */ \
		return x;														\
	}																	\
	SCOPE void tmap_hash_del_##name(tmap_hash_##name##_t *h, tmap_hash_int_t x)				\
	{																	\
		if (x != h->n_buckets && !__tmap_hash_iseither(h->flags, x)) {			\
			__tmap_hash_set_isdel_true(h->flags, x);							\
			--h->size;													\
		}																\
	}

#define TMAP_HASH_INIT(name, tmap_key_t, tmap_val_t, tmap_hash_is_map, __hash_func, __hash_equal) \
	TMAP_HASH_INIT2(name, static inline, tmap_key_t, tmap_val_t, tmap_hash_is_map, __hash_func, __hash_equal)

/* --- BEGIN OF HASH FUNCTIONS --- */

/*!
  Integer hash function
  @param  key   The integer [int32_t]
  @return       The hash value [tmap_hash_int_t]
 */
#define tmap_hash_int_hash_func(key) (int32_t)(key)
/*!
  Integer comparison function
 */
#define tmap_hash_int_hash_equal(a, b) ((a) == (b))
/*!
  64-bit integer hash function
  @param  key   The integer [int64_t]
  @return       The hash value [tmap_hash_int_t]
 */
#define tmap_hash_int64_hash_func(key) (int32_t)((key)>>33^(key)^(key)<<11)
/*!
  64-bit integer comparison function
 */
#define tmap_hash_int64_hash_equal(a, b) ((a) == (b))
/*!
  const char* hash function
  @param  s     Pointer to a null terminated string
  @return       The hash value
 */
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

static inline tmap_hash_int_t __tmap_hash_Wang_hash(tmap_hash_int_t key)
{
    key += ~(key << 15);
    key ^=  (key >> 10);
    key +=  (key << 3);
    key ^=  (key >> 6);
    key += ~(key << 11);
    key ^=  (key >> 16);
    return key;
}
#define tmap_hash_int_hash_func2(k) __tmap_hash_Wang_hash((tmap_hash_int_t)key)

/* --- END OF HASH FUNCTIONS --- */

/* Other convenient macros... */

/*!
  Type of the hash table.
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
  @discussion   For hash sets, calling this results in segfault.
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
#define TMAP_HASH_SET_INIT_INT(name)										\
	TMAP_HASH_INIT(name, int32_t, char, 0, tmap_hash_int_hash_func, tmap_hash_int_hash_equal)

/*!
  Instantiate a hash map containing integer keys
  @param  name  Name of the hash table [symbol]
  @param  tmap_val_t  Type of values [type]
 */
#define TMAP_HASH_MAP_INIT_INT(name, tmap_val_t)								\
	TMAP_HASH_INIT(name, int32_t, tmap_val_t, 1, tmap_hash_int_hash_func, tmap_hash_int_hash_equal)

/*!
  Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
 */
#define TMAP_HASH_SET_INIT_INT64(name)										\
	TMAP_HASH_INIT(name, int64_t, char, 0, tmap_hash_int64_hash_func, tmap_hash_int64_hash_equal)

/*!
  Instantiate a hash map containing 64-bit integer keys
  @param  name  Name of the hash table [symbol]
  @param  tmap_val_t  Type of values [type]
 */
#define TMAP_HASH_MAP_INIT_INT64(name, tmap_val_t)								\
	TMAP_HASH_INIT(name, int64_t, tmap_val_t, 1, tmap_hash_int64_hash_func, tmap_hash_int64_hash_equal)

typedef const char *tmap_hash_cstr_t;
/*!
  Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
 */
#define TMAP_HASH_SET_INIT_STR(name)										\
	TMAP_HASH_INIT(name, tmap_hash_cstr_t, char, 0, tmap_hash_str_hash_func, tmap_hash_str_hash_equal)

/*!
  Instantiate a hash map containing const char* keys
  @param  name  Name of the hash table [symbol]
  @param  tmap_val_t  Type of values [type]
 */
#define TMAP_HASH_MAP_INIT_STR(name, tmap_val_t)								\
	TMAP_HASH_INIT(name, tmap_hash_cstr_t, tmap_val_t, 1, tmap_hash_str_hash_func, tmap_hash_str_hash_equal)

#endif /* __AC_TMAP_HASH_H */
