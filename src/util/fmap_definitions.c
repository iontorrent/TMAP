#include <stdlib.h>
#include <string.h>
#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_definitions.h"
#include "../io/fmap_file.h"

// Input: ASCII character
// Output: 2-bit DNA value
uint8_t fmap_nt_char_to_int[256] = {
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
uint8_t fmap_nt_char_to_rc_char[256] = {
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', '-', 'N', 'N',
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

inline uint32_t fmap_log2(uint32_t v)
{
  uint32_t c = 0;
  if(v & 0xffff0000u) { v >>= 16; c |= 16; }
  if(v & 0xff00) { v >>= 8; c |= 8; }
  if(v & 0xf0) { v >>= 4; c |= 4; }
  if(v & 0xc) { v >>= 2; c |= 2; }
  if(v & 0x2) c |= 1;
  return c;
}

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

  if(0 < end_skip_length) {
      // try first without skipping the end
      char *ptr = fmap_check_suffix(str, suffix, 0);
      if(NULL != ptr) return ptr;
  }

  i=strlen(str) - 1 - end_skip_length;
  j=strlen(suffix)-1;

  if(i < j) return NULL; // the suffix is longer than the string

  while(0 <= j) {
      if(str[i] != suffix[j]) return NULL;
      i--;
      j--;
  }
  return (str + i + 1); // pointer to the start of the suffix
}

static int32_t
fmap_get_last_dot_index(const char *fn)
{
  int32_t len = strlen(fn), l;
  for(l=len-1;0<=l;l--) {
      if('.' == fn[l]) return len-l;
  }
  return len;
}

void
fmap_get_reads_file_format_from_fn_int(char *fn, int32_t *reads_format, int32_t *compr_type)
{
  int32_t compr_suffix_length = 0;

  // auto-recognize the compression type
  if(FMAP_FILE_NO_COMPRESSION == (*compr_type)) {
      if(NULL != fmap_check_suffix(fn, ".bz2", 0)) {
          compr_suffix_length = 4; // ".bz2"
          (*compr_type) = FMAP_FILE_BZ2_COMPRESSION;
      }
      else if(NULL != fmap_check_suffix(fn, ".gz", 0)) {
          compr_suffix_length = 3; // ".gz"
          (*compr_type) = FMAP_FILE_GZ_COMPRESSION;
      }
      else {
          compr_suffix_length = 0; // unknown/none
      }
  }
  else if(FMAP_FILE_BZ2_COMPRESSION == (*compr_type)) {
      compr_suffix_length = 4; // ".bz2"
  }
  else if(FMAP_FILE_GZ_COMPRESSION == (*compr_type)) {
      compr_suffix_length = 3; // ".gz"
  }

  // auto-recognize the reads format
  if(FMAP_READS_FORMAT_UNKNOWN == (*reads_format)) {
      if(NULL != fmap_check_suffix(fn, ".fa", compr_suffix_length) 
         || NULL != fmap_check_suffix(fn, ".fasta", compr_suffix_length)) {
          (*reads_format) = FMAP_READS_FORMAT_FASTA;
      }
      else if(NULL != fmap_check_suffix(fn, ".fq", compr_suffix_length) 
              || NULL != fmap_check_suffix(fn, ".fastq", compr_suffix_length)) {
          (*reads_format) = FMAP_READS_FORMAT_FASTQ;
      }
      else if(NULL != fmap_check_suffix(fn, ".sff", compr_suffix_length)) {
          (*reads_format) = FMAP_READS_FORMAT_SFF;
      }
  }

  // check the suffix implied by the compression type
  switch((*compr_type)) {
    case FMAP_FILE_BZ2_COMPRESSION:
      if(NULL == fmap_check_suffix(fn, ".bz2", 0)) {
          fmap_error("the expected bzip2 file extension is \".bz2\"", Warn, OutOfRange);
          compr_suffix_length = fmap_get_last_dot_index(fn); // remove file extension
      }
      break;
    case FMAP_FILE_GZ_COMPRESSION:
      if(NULL == fmap_check_suffix(fn, ".gz", 0)) {
          fmap_error("the expected gzip file extension is \".gz\"", Warn, OutOfRange);
          compr_suffix_length = fmap_get_last_dot_index(fn); // remove file extension
      }
      break;
    case FMAP_FILE_NO_COMPRESSION:
    default:
      break;
  }

  // check the suffix implied by the reads format
  // Note: try with and without compression file extension
  switch((*reads_format)) {
    case FMAP_READS_FORMAT_FASTA:
      if(NULL == fmap_check_suffix(fn, ".fa", compr_suffix_length) 
         && NULL == fmap_check_suffix(fn, ".fasta", compr_suffix_length)) {
          fmap_error("the expected FASTA file extension is \".fa\" or \".fasta\"", Warn, OutOfRange);
      }
      break;
    case FMAP_READS_FORMAT_FASTQ:
      if(NULL == fmap_check_suffix(fn, ".fq", compr_suffix_length) 
         && NULL == fmap_check_suffix(fn, ".fastq", compr_suffix_length)) {
          fmap_error("the expected FASTA file extension is \".fq\" or \".fastq\"", Warn, OutOfRange);
      }
      break;
    case FMAP_READS_FORMAT_SFF:
      if(NULL == fmap_check_suffix(fn, ".sff", compr_suffix_length)) {
          fmap_error("the expected SFF file extension is \".sff\"", Warn, OutOfRange);
      }
      break;
    case FMAP_READS_FORMAT_UNKNOWN:
    default:
      break;
  }
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

inline void
fmap_reverse(char *seq, int32_t len) 
{
  int32_t i;
  for(i=0;i<(len>>1);i++) {
      char tmp = seq[len-i-1];
      seq[len-i-1] = seq[i];
      seq[i] = tmp;
  }
}

inline void
fmap_reverse_compliment(char *seq, int32_t len) 
{
  int32_t i;
  for(i=0;i<(len>>1);i++) {
      char tmp = seq[len-i-1];
      seq[len-i-1] = fmap_nt_char_to_rc_char[(int)seq[i]];
      seq[i] = fmap_nt_char_to_rc_char[(int)tmp];
  }
  if(1 == (len & 1)) { // mod 2
      seq[i] = fmap_nt_char_to_rc_char[(int)seq[i]];
  }
}
