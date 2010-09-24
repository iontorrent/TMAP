#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_definitions.h"
#include "fmap_string.h"

inline fmap_string_t *
fmap_string_init(int32_t mem)
{
  fmap_string_t *str = NULL;

  str = fmap_calloc(1, sizeof(fmap_string_t), "str");
  if(0 < mem) {
      str->m = mem;
      str->s = fmap_calloc(str->m, sizeof(char), "str->s");
      str->l = 0;
  }

  return str;
}

inline void
fmap_string_destroy(fmap_string_t *str)
{
  if(NULL == str) return;
  free(str->s);
  free(str);
}

inline void
fmap_string_copy(fmap_string_t *dest, fmap_string_t *src)
{
  int32_t i;
  if(dest->m < src->m) {
      dest->m = src->m;
      dest->s = fmap_realloc(dest->s, sizeof(char)*dest->m, "dest->s");
  }
  for(i=0;i<src->m;i++) { // copy over all memory
      dest->s[i] = src->s[i];
  }
  dest->l = src->l;
}

inline fmap_string_t *
fmap_string_clone(fmap_string_t *str)
{
  int32_t i;
  fmap_string_t *ret = NULL;
  
  ret = fmap_string_init(str->m);
  for(i=0;i<str->m;i++) { // copy over all memory
      ret->s[i] = str->s[i];
  }
  ret->l = str->l;

  return ret;
}

inline void
fmap_string_lsprintf(fmap_string_t *dest, int32_t l, const char *format, ...) 
{
  va_list ap;
  int32_t length;
  if(l < 0) fmap_error(NULL, Exit, OutOfRange);
  va_start(ap, format);
  length = vsnprintf(dest->s + l, dest->m - l, format, ap);
  if(length < 0) fmap_error(NULL, Exit, OutOfRange);
  va_end(ap);
  if(dest->m - l - 1 < length) {
      dest->m = length + l + 2;
      fmap_roundup32(dest->m);
      dest->s = fmap_realloc(dest->s, sizeof(char)*dest->m, "dest->s");
      va_start(ap, format);
      length = vsnprintf(dest->s + l, dest->m - l, format, ap);
      va_end(ap);
      if(length < 0) fmap_error(NULL, Exit, OutOfRange);
  }
  dest->l += length;
}

inline void
fmap_string_reverse(fmap_string_t *str)
{
  int i;
  for(i = 0; i < (str->l >> 1); ++i) {
      char tmp = str->s[str->l-1-i];
      str->s[str->l-1-i] = str->s[i]; str->s[i] = tmp;
  }
}

void
fmap_string_reverse_compliment(fmap_string_t *str, int32_t is_int)
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
      for(i = 0; i < (str->l >> 1); ++i) {
          char tmp = str->s[str->l-1-i];
          str->s[str->l-1-i] = nt_char_to_rc_char[(int)str->s[i]]; 
          str->s[i] = nt_char_to_rc_char[(int)tmp];
      }
      if(1 == (str->l & 1)) { // mod 2
          str->s[i] = nt_char_to_rc_char[(int)str->s[i]];
      }
  }
}
