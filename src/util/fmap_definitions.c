#include <stdlib.h>
#include <string.h>
#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_definitions.h"

// Input: ASCII character
// Output: 2-bit DNA value
uint8_t nt_char_to_int[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

inline char *
fmap_get_file_name(const char *prefix, int32_t type)
{
  char *fn = NULL;

  switch(type) {
    case FMAP_ANNO_FILE:
      fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_ANNO_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, FMAP_ANNO_FILE_EXTENSION);
      break;
    case FMAP_PAC_FILE:
      fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_PAC_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, FMAP_PAC_FILE_EXTENSION);
      break;
    case FMAP_REV_PAC_FILE:
      fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_REV_PAC_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, FMAP_REV_PAC_FILE_EXTENSION);
      break;
    case FMAP_BWT_FILE:
      fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_BWT_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, FMAP_BWT_FILE_EXTENSION);
      break;
    case FMAP_REV_BWT_FILE:
      fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_REV_BWT_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, FMAP_REV_BWT_FILE_EXTENSION);
      break;
    case FMAP_SA_FILE:
      fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_SA_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, FMAP_SA_FILE_EXTENSION);
      break;
    case FMAP_REV_SA_FILE:
      fn = fmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(FMAP_REV_SA_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, FMAP_REV_SA_FILE_EXTENSION);
      break;
    default:
      return NULL;
  }
  return fn;
}
