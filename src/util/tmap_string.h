/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_STRING_H_
#define TMAP_STRING_H_

#include <stdint.h>

/*! 
  A Generic String Library
  */

extern uint8_t nt_char_to_rc_char[256];

/*! 
  */
typedef struct {
    size_t l;  /*!< the length of the string */
    size_t m;  /*!< the memory allocated for this string */
    char *s;  /*!< the pointer to the string */
} tmap_string_t;

/*! 
  @param  mem  the initial memory to allocate for the string
  @return      a pointer to the initialized memory
 */
inline tmap_string_t *
tmap_string_init(int32_t mem);

/*! 
  @param  str  a pointer to the string to destroy
 */
inline void
tmap_string_destroy(tmap_string_t *str);

/*! 
  analagous to strcpy
  @param  dest  pointer to the destination string
  @param  src   pointer to the source string
*/
inline void
tmap_string_copy(tmap_string_t *dest, tmap_string_t *src);

/*! 
  @param  str  a pointer to the string to clone
  @return      a pointer to the cloned string
 */
inline tmap_string_t *
tmap_string_clone(tmap_string_t *str);

/*! 
  @param  dest    pointer to the destination string
  @param  l       the number of leading characters to skip
  @param  format  the format for the string
  @param  ...     the arguments to fill in the format
  @details        the first l characters will not be modified
 */
inline void
tmap_string_lsprintf(tmap_string_t *dest, int32_t l, const char *format, ...);

/*! 
  reverse the characters in the string
  @param  str  pointer to the string
*/
inline void
tmap_string_reverse(tmap_string_t *str);

/*! 
  reverse compliments the string
  @param  str       pointer to the string
  @param  is_int    1 if the sequence is in integer format, 0 otherwise
  */
void
tmap_string_reverse_compliment(tmap_string_t *str, int32_t is_int);

#endif
