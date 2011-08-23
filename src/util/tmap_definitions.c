/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <string.h>
#include <config.h>
#include "tmap_error.h"
#include "tmap_alloc.h"
#include "../seq/tmap_seq.h"
#include "../io/tmap_file.h"
#include "tmap_definitions.h"

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
    "mapvsw",
    "mapall"
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

uint8_t tmap_iupac_char_to_int[256] = {
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 0,   5,  1,  6,  15, 15,  2,   7, 15, 15,  8,  15,  9,  4, 15,
    15, 15, 10, 11,  3,  15, 12, 13,  15, 14, 15, 15,  15, 15, 15, 15,
    15, 0,   5,  1,  6,  15, 15,  2,   7, 15, 15,  8,  15,  9,  4, 15,
    15, 15, 10, 11,  3,  15, 12, 13,  15, 14, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15
};

// A -> 1
// C -> 2
// G -> 4
// T -> 8
// R -> 5
// Y -> 10
// S -> 6
// W -> 9
// K -> 12
// M -> 3
// B -> 14
// D -> 13
// H -> 11
// V -> 7
// N -> 15
// A B C D G H K M N R S T V W Y N
uint8_t tmap_iupac_char_to_bit_string[256] = {
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 1,  14,  2,  13, 15, 15,  4,  11, 15, 15, 12,  15,  3, 15, 15,
    15, 15,  5,  6,  8,  15,  7,  9,  15, 10, 15, 15,  15, 15, 15, 15,
    15, 1,  14,  2,  13, 15, 15,  4,  11, 15, 15, 12,  15,  3, 15, 15,
    15, 15,  5,  6,  8,  15,  7,  9,  15, 10, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,
    15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15,  15, 15, 15, 15
};

char tmap_iupac_int_to_char[17] = "ACGTNBDHKMRSVWYN";

int32_t 
tmap_reads_format_to_seq_type(int32_t reads_format)
{
  switch(reads_format) {
    case TMAP_READS_FORMAT_FASTA:
    case TMAP_READS_FORMAT_FASTQ:
      return TMAP_SEQ_TYPE_FQ;
    case TMAP_READS_FORMAT_SFF:
      return TMAP_SEQ_TYPE_SFF;
#ifdef HAVE_SAMTOOLS
    case TMAP_READS_FORMAT_SAM:
      return TMAP_SEQ_TYPE_SAM;
    case TMAP_READS_FORMAT_BAM:
      return TMAP_SEQ_TYPE_BAM;
#endif
    default:
      return TMAP_SEQ_TYPE_NOTYPE;
  }
}


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
#ifdef HAVE_SAMTOOLS
  else if(0 == strcmp(optarg, "sam")) {
      return TMAP_READS_FORMAT_SAM;
  }
  else if(0 == strcmp(optarg, "bam")) {
      return TMAP_READS_FORMAT_BAM;
  }
#endif
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
      if(NULL != tmap_check_suffix(fn, ".gz", 0)) {
          compr_suffix_length = 3; // ".gz"
          (*compr_type) = TMAP_FILE_GZ_COMPRESSION;
      }
#ifndef DISABLE_BZ2
      else if(NULL != tmap_check_suffix(fn, ".bz2", 0)) {
          compr_suffix_length = 4; // ".bz2"
          (*compr_type) = TMAP_FILE_BZ2_COMPRESSION;
      }
#endif
      else {
          compr_suffix_length = 0; // unknown/none
      }
  }
#ifndef DISABLE_BZ2
  else if(TMAP_FILE_BZ2_COMPRESSION == (*compr_type)) {
      compr_suffix_length = 4; // ".bz2"
  }
#endif
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
#ifdef HAVE_SAMTOOLS
      else if(NULL != tmap_check_suffix(fn, ".sam", compr_suffix_length)) {
          (*reads_format) = TMAP_READS_FORMAT_SAM;
      }
      else if(NULL != tmap_check_suffix(fn, ".bam", compr_suffix_length)) {
          (*reads_format) = TMAP_READS_FORMAT_BAM;
      }
#endif
  }

  // check the suffix implied by the compression type
  switch((*compr_type)) {
#ifndef DISABLE_BZ2
    case TMAP_FILE_BZ2_COMPRESSION:
      if(NULL == tmap_check_suffix(fn, ".bz2", 0)) {
          tmap_error("the expected bzip2 file extension is \".bz2\"", Warn, OutOfRange);
          compr_suffix_length = tmap_get_last_dot_index(fn); // remove file extension
      }
      break;
#endif
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
#ifdef HAVE_SAMTOOLS
    case TMAP_READS_FORMAT_SAM:
      if(NULL == tmap_check_suffix(fn, ".sam", compr_suffix_length)) {
          tmap_error("the expected SAM file extension is \".sam\"", Warn, OutOfRange);
      }
      break;
    case TMAP_READS_FORMAT_BAM:
      if(NULL == tmap_check_suffix(fn, ".bam", compr_suffix_length)) {
          tmap_error("the expected BAM file extension is \".bam\"", Warn, OutOfRange);
      }
      break;
#endif
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
#ifdef HAVE_SAMTOOLS
  else if(TMAP_READS_FORMAT_SAM == format) {
      return strdup("sam");
  }
  else if(TMAP_READS_FORMAT_BAM == format) {
      return strdup("bam");
  }
#endif
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
tmap_reverse_int(uint8_t *seq, int32_t len) 
{
  int32_t i;
  for(i=0;i<(len>>1);i++) {
      uint8_t tmp = seq[len-i-1];
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
tmap_reverse_compliment_int(uint8_t *seq, int32_t len) 
{
  int32_t i;
  for(i=0;i<(len>>1);i++) {
      char tmp = seq[len-i-1];
      seq[len-i-1] = (4 <= seq[i]) ? seq[i] : (3 - seq[i]);
      seq[i] = (4 <= tmp) ? tmp : (3 - tmp);
  }
  if(1 == (len & 1)) { // mod 2
      seq[i] = (4 <= seq[i]) ? seq[i] : (3 - seq[i]);
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

int32_t
tmap_chomp(char *str)
{
  int32_t n, len;

  if(NULL == str) return 0;

  len = strlen(str);
  for(n=len;0<n;n--) {
      if(0 == isspace(str[n-1])) break;
      str[n-1]='\0';
  }

  return len - n;
}

inline int32_t
tmap_interval_overlap(uint32_t low1, uint32_t high1, uint32_t low2, uint32_t high2)
{
  if(high1 < low2) {
      return -1;
  }
  else if(high2 < low1) {
      return 1;
  }
  else {
      return 0;
  }
}


static inline void
tmap_version_to_int(const char *v, int32_t v_n[3])
{
  int32_t i, j, len;
  len = strlen(v);
  v_n[0] = atoi(v);
  for(i=0,j=1;i<len;i++) {
      if('.' == v[i]) {
          i++; // skip the dot
          if(i == len) {
              tmap_error("malformed version string", Exit, OutOfRange);
          }
          v_n[j] = atoi(v + i);
          j++;
      }
  }
  if(3 != j) {
      tmap_error("malformed version string", Exit, OutOfRange);
  }
}

int32_t
tmap_compare_versions(const char *v1, const char *v2)
{
  int32_t i;
  int32_t v1_n[3], v2_n[3];

  tmap_version_to_int(v1, v1_n);
  tmap_version_to_int(v2, v2_n);
  
  for(i=0;i<3;i++) {
      if(v1_n[i] < v2_n[i]) return -1;
      else if(v1_n[i] > v2_n[i]) return 1;
  }

  return 0;
}

int32_t
tmap_validate_flow_order(const char *flow_order)
{
  int32_t i, j;
  tmap_error_cmd_check_int(strlen(flow_order), 4, INT32_MAX, "-x");
  for(i=j=0;i<strlen(flow_order);i++) { // each base must be used
      switch(tolower(flow_order[i])) {
        case 'a': j |= 0x1; break;
        case 'c': j |= 0x2; break;
        case 'g': j |= 0x4; break;
        case 't': j |= 0x8; break;
        default: return -1; 
      }
  }
  if(0xf != j) {
      return -2;
  }
  return 0;
}

int32_t
tmap_validate_key_seq(const char *key_seq)
{
  int32_t i;
  tmap_error_cmd_check_int(strlen(key_seq), 4, INT32_MAX, "-x");
  for(i=0;i<strlen(key_seq);i++) { // each base must be used
      switch(tolower(key_seq[i])) {
        case 'a': 
        case 'c': 
        case 'g': 
        case 't': 
          break;
        default: 
          return -1;
      }
  }
  return 0;
}
