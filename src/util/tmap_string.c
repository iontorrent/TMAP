/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include "tmap_error.h"
#include "tmap_alloc.h"
#include "tmap_definitions.h"
#include "tmap_string.h"

inline tmap_string_t *
tmap_string_init(int32_t mem)
{
  tmap_string_t *str = NULL;

  str = tmap_calloc(1, sizeof(tmap_string_t), "str");
  if(0 < mem) {
      str->m = mem;
      str->s = tmap_calloc(str->m, sizeof(char), "str->s");
      str->l = 0;
      str->s[str->l] = '\0';
  }

  return str;
}

inline void
tmap_string_destroy(tmap_string_t *str)
{
  if(NULL == str) return;
  free(str->s);
  free(str);
}

inline void
tmap_string_copy(tmap_string_t *dest, tmap_string_t *src)
{
  int32_t i;
  if(dest->m < src->m) {
      dest->m = src->m;
      dest->s = tmap_realloc(dest->s, sizeof(char)*dest->m, "dest->s");
  }
  for(i=0;i<src->m;i++) { // copy over all memory
      dest->s[i] = src->s[i];
  }
  dest->l = src->l;
}

inline tmap_string_t *
tmap_string_clone(tmap_string_t *str)
{
  int32_t i;
  tmap_string_t *ret = NULL;
  
  ret = tmap_string_init(str->m);
  for(i=0;i<str->m;i++) { // copy over all memory
      ret->s[i] = str->s[i];
  }
  ret->l = str->l;

  return ret;
}

inline void
tmap_string_lsprintf(tmap_string_t *dest, int32_t l, const char *format, ...) 
{
  va_list ap;
  int32_t length;
  if(l < 0) tmap_error(NULL, Exit, OutOfRange);
  va_start(ap, format);
  length = vsnprintf(dest->s + l, dest->m - l, format, ap);
  if(length < 0) tmap_error(NULL, Exit, OutOfRange);
  va_end(ap);
  if(dest->m - l - 1 < length) {
      dest->m = length + l + 2;
      tmap_roundup32(dest->m);
      dest->s = tmap_realloc(dest->s, sizeof(char)*dest->m, "dest->s");
      va_start(ap, format);
      length = vsnprintf(dest->s + l, dest->m - l, format, ap);
      va_end(ap);
      if(length < 0) tmap_error(NULL, Exit, OutOfRange);
  }
  dest->l += length;
}

inline void
tmap_string_reverse(tmap_string_t *str)
{
  tmap_reverse(str->s, str->l);
}

void
tmap_string_reverse_compliment(tmap_string_t *str, int32_t is_int)
{
  int i;

  if(1 == is_int) { // bases are integer values
      for(i = 0; i < (str->l >> 1); ++i) {
          char tmp = str->s[str->l-1-i];
          str->s[str->l-1-i] = (4 <= str->s[i]) ? str->s[i] : 3 - str->s[i];
          str->s[i] = (4 <= tmp) ? tmp : 3 - tmp;
      }
      if(1 == (str->l & 1)) { // mod 2
          str->s[i] = (4 <= str->s[i]) ? str->s[i] : 3 - str->s[i];
      }
  }
  else { // bases are ASCII values
      tmap_reverse_compliment(str->s, str->l);
  }
}
