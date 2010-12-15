/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_ALLOC_H
#define TMAP_ALLOC_H

#include "tmap_error.h"

/*! 
  Memory Allocation Routines. 
  */

/*! 
  wrapper function for malloc
  @param  _size           the size of the memory block, in bytes
  @param  _variable_name  the variable name to be assigned this memory in the calling function
  @return                 upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  */
#define tmap_malloc(_size, _variable_name) \
  tmap_malloc1(_size, __func__, _variable_name)

/*! 
  wrapper function for realloc
  @param  _ptr            the variable to be reallocated
  @param  _size           the size of the memory block, in bytes
  @param  _variable_name  the variable name to be assigned this memory in the calling function
  @return 		  upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  @details		  the ptr must be a memory block previously allocated with malloc, calloc, or realloc to be reallocated; if the ptr is NULL, a new block of memory will be allocated. 
  */
#define tmap_realloc(_ptr, _size, _variable_name) \
  tmap_realloc1(_ptr, _size, __func__, _variable_name)

/*! 
  wrapper function for calloc
  @param  _num            the number of elements to be allocated
  @param  _size           the size of the memory block, in bytes
  @param  _variable_name  the variable name to be assigned this memory in the calling function
  @return                 upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  */
#define tmap_calloc(_num, _size, _variable_name) \
  tmap_calloc1(_num, _size, __func__, _variable_name)

/*! 
  wrapper for 'strdup' that checks memory allocation
  @param  _str            string to be copied
  @return                 a pointer to the copied string
  */
#define tmap_strdup(_str) \
  tmap_strdup1(_str, __func__)

/*! 
  wrapper function for malloc
  @param  size           the size of the memory block, in bytes
  @param  function_name  the calling function name 
  @param  variable_name  the variable name to be assigned this memory in the calling function
  @return                upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  */
inline void *
tmap_malloc1(size_t size, const char *function_name, const char *variable_name);

/*! 
  wrapper function for realloc
  @param  ptr            the variable to be reallocated
  @param  size           the size of the memory block, in bytes
  @param  function_name  the calling function name 
  @param  variable_name  the variable name to be assigned this memory in the calling function
  @return 		 upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  @details		 the ptr must be a memory block previously allocated with malloc, calloc, or realloc to be reallocated; if the ptr is NULL, a new block of memory will be allocated. 
  */
inline void *
tmap_realloc1(void *ptr, size_t size, const char *function_name, const char *variable_name);

/*! 
  wrapper function for calloc
  @param  num            the number of elements to be allocated
  @param  size           the size of the memory block, in bytes
  @param  function_name  the calling function name 
  @param  variable_name  the variable name to be assigned this memory in the calling function
  @return                upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  */
inline void *
tmap_calloc1(size_t num, size_t size, const char *function_name, const char *variable_name);

/*! 
  wrapper for 'strdup' that checks memory allocation
  @param  str            string to be copied
  @param  function_name  the calling function name 
  @return                a pointer to the copied string
  */
inline char *
tmap_strdup1(const char *str, const char *function_name);

#endif
