#include <stdlib.h>
#include <string.h>
#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_definitions.h"
#include "../io/fmap_file.h"

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

// Input: ASCII character
// Output: ASCII reverse complimented DNA value
uint8_t nt_char_to_rc_char[256] = {
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'A', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G',  'N', 'N', 'N', 'C',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'A', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
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

int
fmap_get_reads_file_format_int(char *optarg)
{
  if(0 == strcmp(optarg, "fa") || 0 == strcmp(optarg, "fasta")) {
      return FMAP_READS_FORMAT_FASTA;
  }
  else if(0 == strcmp(optarg, "fq") || 0 == strcmp(optarg, "fastq")) {
      return FMAP_READS_FORMAT_FASTQ;
  }
  else if(0 == strcmp(optarg, "sff")) {
      return FMAP_READS_FORMAT_SFF;
  }
  return FMAP_READS_FORMAT_UNKNOWN;
}

static char *
fmap_check_suffix(char *str, char *suffix, int32_t end_skip_length)
{
  int32_t i, j;

  if(NULL == str) return NULL;
  if(NULL == suffix) return str; // empty suffix always matches

  i=strlen(str)- 1 - end_skip_length;
  j=strlen(suffix)-1;

  if(i < j) return NULL; // the suffix is longer than the string

  while(0 <= j) {
      if(str[i] != suffix[j]) return NULL;
      i--;
      j--;
  }
  return (str + i + 1); // pointer to the start of the suffix
}

int
fmap_get_reads_file_format_from_fn_int(char *fn, int32_t compr_type)
{
  int32_t compr_suffix_length = 0;
  switch(compr_type) {
    case FMAP_FILE_BZ2_COMPRESSION:
      if(NULL == fmap_check_suffix(fn, ".bz2", 0)) {
          fmap_error("the expected bzip2 file extension is \".bz2\"", Warn, OutOfRange);
          return FMAP_READS_FORMAT_UNKNOWN;
      }
      compr_suffix_length = 4; // ".bz2"
      break;
    case FMAP_FILE_GZ_COMPRESSION:
      if(NULL == fmap_check_suffix(fn, ".gz", 0)) {
          fmap_error("the expected gzip file extension is \".gz\"", Warn, OutOfRange);
          return FMAP_READS_FORMAT_UNKNOWN;
      }
      compr_suffix_length = 3; // ".gz"
      break;
    case FMAP_FILE_NO_COMPRESSION:
    default:
      break;
  }
  if(NULL != fmap_check_suffix(fn, ".fa", compr_suffix_length) 
     || NULL != fmap_check_suffix(fn, ".fasta", compr_suffix_length)) {
      return FMAP_READS_FORMAT_FASTA;
  }
  else if(NULL != fmap_check_suffix(fn, ".fq", compr_suffix_length) 
          || NULL != fmap_check_suffix(fn, ".fastq", compr_suffix_length)) {
      return FMAP_READS_FORMAT_FASTQ;
  }
  else if(NULL != fmap_check_suffix(fn, ".sff", compr_suffix_length)) {
      return FMAP_READS_FORMAT_SFF;
  }
  return FMAP_READS_FORMAT_UNKNOWN;
}

char *
fmap_get_reads_file_format_string(int format)
{
  if(FMAP_READS_FORMAT_FASTA == format) {
      return strdup("fasta");
  }
  else if(FMAP_READS_FORMAT_FASTQ == format) {
      return strdup("fastq");
  }
  else if(FMAP_READS_FORMAT_SFF == format) {
      return strdup("sff");
  }
  else if(FMAP_READS_FORMAT_UNKNOWN == format) {
      return strdup("unknown");
  }
  else {
      return strdup("unknown1");
  }

}
