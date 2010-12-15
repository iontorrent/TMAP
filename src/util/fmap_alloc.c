#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "fmap_alloc.h"

inline void *
fmap_malloc1(size_t size, const char *function_name, const char *variable_name)
{
  void *ptr = malloc(size);
  if(NULL == ptr) {
      fmap_error1(function_name, variable_name, Exit, MallocMemory);
  }
  return ptr;
}

inline void *
fmap_realloc1(void *ptr, size_t size, const char *function_name, const char *variable_name)
{
  ptr = realloc(ptr, size);
  if(0 != size && NULL == ptr) {
      fmap_error1(function_name, variable_name, Exit, ReallocMemory);
  }
  return ptr;
}

inline void *
fmap_calloc1(size_t num, size_t size, const char *function_name, const char *variable_name)
{
  void *ptr = calloc(num, size);
  if(NULL == ptr) {
      fmap_error1(function_name, variable_name, Exit, MallocMemory);
  }
  return ptr;
}

inline char *
fmap_strdup1(const char *str, const char *function_name)
{
  if(NULL == str) return NULL;
  char *ptr = strdup(str);
  if(NULL == ptr) {
      fmap_error1(function_name, NULL, Exit, MallocMemory);
  }
  return ptr;
}
