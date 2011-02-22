/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <string.h>
#include "tmap_error.h"
#include "tmap_alloc.h"
#include "tmap_definitions.h"
#include "../io/tmap_file.h"

// Algorithm IDs

static char *algo_id_to_name[17] = {
    "none",
    "map1", 
    "map2", 
    "map3", 
    "dummy4",
    "dummy5",
    "dummy6",
    "dummy7",
    "dummy8",
    "dummy9",
    "dummy10",
    "dummy11",
    "dummy12",
    "dummy13",
    "dummy14",
    "mapall",
    "mappability"
};

char *
tmap_algo_id_to_name(uint16_t algo_id)
{
  int32_t i=0;
  while(0 < algo_id) {
      algo_id >>= 1;
      i++;
  }
  return algo_id_to_name[i];
}

// Input: ASCII character
// Output: 2-bit DNA value
uint8_t tmap_nt_char_to_int[256] = {
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
uint8_t tmap_nt_char_to_rc_char[256] = {
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

inline uint32_t tmap_log2(uint32_t v)
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
tmap_get_file_name(const char *prefix, int32_t type)
{
  char *fn = NULL;

  switch(type) {
    case TMAP_ANNO_FILE:
      fn = tmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(TMAP_ANNO_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, TMAP_ANNO_FILE_EXTENSION);
      break;
    case TMAP_PAC_FILE:
      fn = tmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(TMAP_PAC_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, TMAP_PAC_FILE_EXTENSION);
      break;
    case TMAP_REV_PAC_FILE:
      fn = tmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(TMAP_REV_PAC_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, TMAP_REV_PAC_FILE_EXTENSION);
      break;
    case TMAP_BWT_FILE:
      fn = tmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(TMAP_BWT_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, TMAP_BWT_FILE_EXTENSION);
      break;
    case TMAP_REV_BWT_FILE:
      fn = tmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(TMAP_REV_BWT_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, TMAP_REV_BWT_FILE_EXTENSION);
      break;
    case TMAP_SA_FILE:
      fn = tmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(TMAP_SA_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, TMAP_SA_FILE_EXTENSION);
      break;
    case TMAP_REV_SA_FILE:
      fn = tmap_malloc(sizeof(char)*(1+strlen(prefix)+strlen(TMAP_REV_SA_FILE_EXTENSION)), "fn");
      strcpy(fn, prefix);
      strcat(fn, TMAP_REV_SA_FILE_EXTENSION);
      break;
    default:
      return NULL;
  }
  return fn;
}

int
tmap_get_reads_file_format_int(char *optarg)
{
  if(0 == strcmp(optarg, "fa") || 0 == strcmp(optarg, "fasta")) {
      return TMAP_READS_FORMAT_FASTA;
  }
  else if(0 == strcmp(optarg, "fq") || 0 == strcmp(optarg, "fastq")) {
      return TMAP_READS_FORMAT_FASTQ;
  }
  else if(0 == strcmp(optarg, "sff")) {
      return TMAP_READS_FORMAT_SFF;
  }
  return TMAP_READS_FORMAT_UNKNOWN;
}

static char *
tmap_check_suffix(char *str, char *suffix, int32_t end_skip_length)
{
  int32_t i, j;

  if(NULL == str) return NULL;
  if(NULL == suffix) return str; // empty suffix always matches

  if(0 < end_skip_length) {
      // try first without skipping the end
      char *ptr = tmap_check_suffix(str, suffix, 0);
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
tmap_get_last_dot_index(const char *fn)
{
  int32_t len = strlen(fn), l;
  for(l=len-1;0<=l;l--) {
      if('.' == fn[l]) return len-l;
  }
  return len;
}

void
tmap_get_reads_file_format_from_fn_int(char *fn, int32_t *reads_format, int32_t *compr_type)
{
  int32_t compr_suffix_length = 0;

  // auto-recognize the compression type
  if(TMAP_FILE_NO_COMPRESSION == (*compr_type)) {
      if(NULL != tmap_check_suffix(fn, ".bz2", 0)) {
          compr_suffix_length = 4; // ".bz2"
          (*compr_type) = TMAP_FILE_BZ2_COMPRESSION;
      }
      else if(NULL != tmap_check_suffix(fn, ".gz", 0)) {
          compr_suffix_length = 3; // ".gz"
          (*compr_type) = TMAP_FILE_GZ_COMPRESSION;
      }
      else {
          compr_suffix_length = 0; // unknown/none
      }
  }
  else if(TMAP_FILE_BZ2_COMPRESSION == (*compr_type)) {
      compr_suffix_length = 4; // ".bz2"
  }
  else if(TMAP_FILE_GZ_COMPRESSION == (*compr_type)) {
      compr_suffix_length = 3; // ".gz"
  }

  if(NULL == fn) return;

  // auto-recognize the reads format
  if(TMAP_READS_FORMAT_UNKNOWN == (*reads_format)) {
      if(NULL != tmap_check_suffix(fn, ".fa", compr_suffix_length) 
         || NULL != tmap_check_suffix(fn, ".fasta", compr_suffix_length)) {
          (*reads_format) = TMAP_READS_FORMAT_FASTA;
      }
      else if(NULL != tmap_check_suffix(fn, ".fq", compr_suffix_length) 
              || NULL != tmap_check_suffix(fn, ".fastq", compr_suffix_length)) {
          (*reads_format) = TMAP_READS_FORMAT_FASTQ;
      }
      else if(NULL != tmap_check_suffix(fn, ".sff", compr_suffix_length)) {
          (*reads_format) = TMAP_READS_FORMAT_SFF;
      }
  }

  // check the suffix implied by the compression type
  switch((*compr_type)) {
    case TMAP_FILE_BZ2_COMPRESSION:
      if(NULL == tmap_check_suffix(fn, ".bz2", 0)) {
          tmap_error("the expected bzip2 file extension is \".bz2\"", Warn, OutOfRange);
          compr_suffix_length = tmap_get_last_dot_index(fn); // remove file extension
      }
      break;
    case TMAP_FILE_GZ_COMPRESSION:
      if(NULL == tmap_check_suffix(fn, ".gz", 0)) {
          tmap_error("the expected gzip file extension is \".gz\"", Warn, OutOfRange);
          compr_suffix_length = tmap_get_last_dot_index(fn); // remove file extension
      }
      break;
    case TMAP_FILE_NO_COMPRESSION:
    default:
      break;
  }

  // check the suffix implied by the reads format
  // Note: try with and without compression file extension
  switch((*reads_format)) {
    case TMAP_READS_FORMAT_FASTA:
      if(NULL == tmap_check_suffix(fn, ".fa", compr_suffix_length) 
         && NULL == tmap_check_suffix(fn, ".fasta", compr_suffix_length)) {
          tmap_error("the expected FASTA file extension is \".fa\" or \".fasta\"", Warn, OutOfRange);
      }
      break;
    case TMAP_READS_FORMAT_FASTQ:
      if(NULL == tmap_check_suffix(fn, ".fq", compr_suffix_length) 
         && NULL == tmap_check_suffix(fn, ".fastq", compr_suffix_length)) {
          tmap_error("the expected FASTA file extension is \".fq\" or \".fastq\"", Warn, OutOfRange);
      }
      break;
    case TMAP_READS_FORMAT_SFF:
      if(NULL == tmap_check_suffix(fn, ".sff", compr_suffix_length)) {
          tmap_error("the expected SFF file extension is \".sff\"", Warn, OutOfRange);
      }
      break;
    case TMAP_READS_FORMAT_UNKNOWN:
    default:
      break;
  }
}

char *
tmap_get_reads_file_format_string(int format)
{
  if(TMAP_READS_FORMAT_FASTA == format) {
      return strdup("fasta");
  }
  else if(TMAP_READS_FORMAT_FASTQ == format) {
      return strdup("fastq");
  }
  else if(TMAP_READS_FORMAT_SFF == format) {
      return strdup("sff");
  }
  else if(TMAP_READS_FORMAT_UNKNOWN == format) {
      return strdup("unknown");
  }
  else {
      return strdup("unknown1");
  }

}

inline void
tmap_reverse(char *seq, int32_t len) 
{
  int32_t i;
  for(i=0;i<(len>>1);i++) {
      char tmp = seq[len-i-1];
      seq[len-i-1] = seq[i];
      seq[i] = tmp;
  }
}

inline void
tmap_reverse_compliment(char *seq, int32_t len) 
{
  int32_t i;
  for(i=0;i<(len>>1);i++) {
      char tmp = seq[len-i-1];
      seq[len-i-1] = tmap_nt_char_to_rc_char[(int)seq[i]];
      seq[i] = tmap_nt_char_to_rc_char[(int)tmp];
  }
  if(1 == (len & 1)) { // mod 2
      seq[i] = tmap_nt_char_to_rc_char[(int)seq[i]];
  }
}

inline void
tmap_compliment(char *seq, int32_t len) 
{
  int32_t i;
  for(i=0;i<len;i++) {
      seq[i] = tmap_nt_char_to_rc_char[(int)seq[i]];
  }
}
