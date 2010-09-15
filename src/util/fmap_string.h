#ifndef FMAP_STRING_H_
#define FMAP_STRING_H_

#include <stdint.h>

extern uint8_t nt_char_to_rc_char[256];

/*! @typedef
  @abstract
  @field  l  the length of the string
  @field  m  the memory allocated for this string
  @field  s  the pointer to the string
  */
typedef struct {
    size_t l, m;
    char *s;
} fmap_string_t;

/*! @function
  @param  mem  the initial memory to allocate for the string
  @return      a pointer to the initialized memory
 */
inline fmap_string_t *
fmap_string_init(int32_t mem);

/*! @function
  @abstract
  @param  str  a pointer to the string to destroy
 */
inline void
fmap_string_destroy(fmap_string_t *str);

/*! @function
  @abstract     analagous to strcpy
  @param  dest  pointer to the destination string
  @param  src   pointer to the source string
*/
inline void
fmap_string_copy(fmap_string_t *dest, fmap_string_t *src);

/*! @function
  @abstract
  @param  str  a pointer to the string to clone
  @return      a pointer to the cloned string
 */
inline fmap_string_t *
fmap_string_clone(fmap_string_t *str);

/*! @function
  @abstract    reverse the characters in the string
  @param  str  pointer to the string
*/
inline void
fmap_string_reverse(fmap_string_t *str);

/*! @function
  @abstract         reverse compliments the string
  @param  str       pointer to the string
  @param  is_int    1 if the sequence is in integer format, 0 otherwise
  */
void
fmap_string_reverse_compliment(fmap_string_t *str, int32_t is_int);

#endif
