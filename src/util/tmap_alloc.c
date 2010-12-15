/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "tmap_alloc.h"

inline void *
tmap_malloc1(size_t size, const char *function_name, const char *variable_name)
{
  void *ptr = malloc(size);
  if(NULL == ptr) {
      tmap_error1(function_name, variable_name, Exit, MallocMemory);
  }
  return ptr;
}

inline void *
tmap_realloc1(void *ptr, size_t size, const char *function_name, const char *variable_name)
{
  ptr = realloc(ptr, size);
  if(0 != size && NULL == ptr) {
      tmap_error1(function_name, variable_name, Exit, ReallocMemory);
  }
  return ptr;
}

inline void *
tmap_calloc1(size_t num, size_t size, const char *function_name, const char *variable_name)
{
  void *ptr = calloc(num, size);
  if(NULL == ptr) {
      tmap_error1(function_name, variable_name, Exit, MallocMemory);
  }
  return ptr;
}

inline char *
tmap_strdup1(const char *str, const char *function_name)
{
  if(NULL == str) return NULL;
  char *ptr = strdup(str);
  if(NULL == ptr) {
      tmap_error1(function_name, NULL, Exit, MallocMemory);
  }
  return ptr;
}
