#ifndef FMAP_ALLOC_H_
#define FMAP_ALLOC_H_

#include "fmap_error.h"

/*! 
  Memory Allocation Routines. 
  */

/*! 
               wrapper function for malloc
  @param  _size           the size of the memory block, in bytes
  @param  _variable_name  the variable name to be assigned this memory in the calling function
  @return                 upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  */
#define fmap_malloc(_size, _variable_name) \
  fmap_malloc1(_size, __func__, _variable_name)

/*! 
               wrapper function for realloc
  @param  _ptr            the variable to be reallocated
  @param  _size           the size of the memory block, in bytes
  @param  _variable_name  the variable name to be assigned this memory in the calling function
  @return 		  upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  details		  the ptr must be a memory block previously allocated with malloc, calloc, or realloc to be reallocated; if the ptr is NULL, a new block of memory will be allocated. 
  */
#define fmap_realloc(_ptr, _size, _variable_name) \
  fmap_realloc1(_ptr, _size, __func__, _variable_name)

/*! 
               wrapper function for calloc
  @param  _num            the number of elements to be allocated
  @param  _size           the size of the memory block, in bytes
  @param  _variable_name  the variable name to be assigned this memory in the calling function
  @return                 upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  */
#define fmap_calloc(_num, _size, _variable_name) \
  fmap_calloc1(_num, _size, __func__, _variable_name)

/*! 
               wrapper for 'strdup' that checks memory allocation
  @param  _str            string to be copied
  @return                 a pointer to the copied string
  */
#define fmap_strdup(_str) \
  fmap_strdup1(_str, __func__)

/*! 
              wrapper function for malloc
  @param  size           the size of the memory block, in bytes
  @param  function_name  the calling function name 
  @param  variable_name  the variable name to be assigned this memory in the calling function
  @return                upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  */
inline void *
fmap_malloc1(size_t size, const char *function_name, const char *variable_name);

/*! 
              wrapper function for realloc
  @param  ptr            the variable to be reallocated
  @param  size           the size of the memory block, in bytes
  @param  function_name  the calling function name 
  @param  variable_name  the variable name to be assigned this memory in the calling function
  @return 		 upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  details		 the ptr must be a memory block previously allocated with malloc, calloc, or realloc to be reallocated; if the ptr is NULL, a new block of memory will be allocated. 
  */
inline void *
fmap_realloc1(void *ptr, size_t size, const char *function_name, const char *variable_name);

/*! 
  wrapper function for calloc
  @param  num            the number of elements to be allocated
  @param  size           the size of the memory block, in bytes
  @param  function_name  the calling function name 
  @param  variable_name  the variable name to be assigned this memory in the calling function
  @return                upon success, a pointer to the memory block allocated by the function; a null pointer otherwise.
  */
inline void *
fmap_calloc1(size_t num, size_t size, const char *function_name, const char *variable_name);

/*! 
              wrapper for 'strdup' that checks memory allocation
  @param  str            string to be copied
  @param  function_name  the calling function name 
  @return                a pointer to the copied string
  */
inline char *
fmap_strdup1(const char *str, const char *function_name);

#endif
