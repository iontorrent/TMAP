/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include <config.h>
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_error.h"
#include "../../io/tmap_file.h"
#include "../../seq/tmap_seq.h"
#include "../../util/tmap_progress.h"
#include "../../util/tmap_definitions.h"
#include "tmap_map_opt.h"
  
static char *tmap_map_opt_input_types[] = {"INT", "FLOAT", "NUM", "FILE", "STRING", "NONE"};

// int32_t print function
#define __tmap_map_opt_option_print_func_int_init(_name) \
  static void tmap_map_opt_option_print_func_##_name(void *arg) { \
      tmap_map_opt_t *opt = (tmap_map_opt_t*)arg; \
      tmap_file_fprintf(tmap_file_stderr, "[%d]", opt->_name); \
  }

// double print function
#define __tmap_map_opt_option_print_func_double_init(_name) \
  static void tmap_map_opt_option_print_func_##_name(void *arg) { \
      tmap_map_opt_t *opt = (tmap_map_opt_t*)arg; \
      tmap_file_fprintf(tmap_file_stderr, "[%lf]", opt->_name); \
  }

// char array print function
#define __tmap_map_opt_option_print_func_chars_init(_name, _null_msg) \
  static void tmap_map_opt_option_print_func_##_name(void *arg) { \
      tmap_map_opt_t *opt = (tmap_map_opt_t*)arg; \
      tmap_file_fprintf(tmap_file_stderr, "[%s]", (NULL == opt->_name) ? _null_msg : opt->_name); \
  } \

// int32_t print function
#define __tmap_map_opt_option_print_func_int_autodetected_init(_name, _detected) \
  static void tmap_map_opt_option_print_func_##_name(void *arg) { \
      tmap_map_opt_t *opt = (tmap_map_opt_t*)arg; \
      tmap_file_fprintf(tmap_file_stderr, "(%s) [%d]", (0 == opt->_detected) ? "user set" : "autodetect", opt->_name); \
  }

// true/false 
#define __tmap_map_opt_option_print_func_tf_init(_name) \
  static void tmap_map_opt_option_print_func_##_name(void *arg) { \
      tmap_map_opt_t *opt = (tmap_map_opt_t*)arg; \
      tmap_file_fprintf(tmap_file_stderr, "[%s]", (1 == opt->_name) ? "true" : "false"); \
  }

// verbosity
#define __tmap_map_opt_option_print_func_verbosity_init() \
  static void tmap_map_opt_option_print_func_verbosity(void *arg) { \
      tmap_file_fprintf(tmap_file_stderr, "[%s]", (1 == tmap_progress_get_verbosity()) ? "true" : "false"); \
  }

// compression
#define __tmap_map_opt_option_print_func_compr_init(_name, _var, _compr) \
  static void tmap_map_opt_option_print_func_##_name(void *arg) { \
      tmap_map_opt_t *opt = (tmap_map_opt_t*)arg; \
      tmap_file_fprintf(tmap_file_stderr, "[%s]", (_compr == opt->_var) ? "using" : "not using"); \
  }

// number/probability
#define __tmap_map_opt_option_print_func_np_init(_name_num, _name_prob) \
  static void tmap_map_opt_option_print_func_##_name_num(void *arg) { \
      tmap_map_opt_t *opt = (tmap_map_opt_t*)arg; \
      if(opt->_name_num < 0) { \
          tmap_file_fprintf(tmap_file_stderr, "[probability: %lf]", opt->_name_prob); \
      } \
      else { \
          tmap_file_fprintf(tmap_file_stderr, "[number: %d]", opt->_name_num); \
      } \
  }

// number/fraction
#define __tmap_map_opt_option_print_func_nf_init(_name_num, _name_frac) \
  static void tmap_map_opt_option_print_func_##_name_num(void *arg) { \
      tmap_map_opt_t *opt = (tmap_map_opt_t*)arg; \
      if(opt->_name_num < 0) { \
          tmap_file_fprintf(tmap_file_stderr, "[fraction: %lf]", opt->_name_frac); \
      } \
      else { \
          tmap_file_fprintf(tmap_file_stderr, "[number: %d]", opt->_name_num); \
      } \
  }

// array of character arrays (i.e. list of file names)
#define __tmap_map_opt_option_print_func_char_array_init(_name_array, _name_length, _default) \
  static void tmap_map_opt_option_print_func_##_name_array(void *arg) { \
      tmap_map_opt_t *opt = (tmap_map_opt_t*)arg; \
      int32_t i; \
      tmap_file_fprintf(tmap_file_stderr, "["); \
      if(0 == opt->_name_length) tmap_file_fprintf(tmap_file_stderr, _default); \
      else { \
          for(i=0;i<opt->_name_length;i++) { \
              if(0 < i) tmap_file_fprintf(tmap_file_stderr, ","); \
              tmap_file_fprintf(tmap_file_stderr, "%s", opt->_name_array[i]); \
          } \
      } \
      tmap_file_fprintf(tmap_file_stderr, "]"); \
  }

// reads format
#define __tmap_map_opt_option_print_func_reads_format_init(_name) \
  static void tmap_map_opt_option_print_func_##_name(void *arg) { \
      tmap_map_opt_t *opt = (tmap_map_opt_t*)arg; \
      char *reads_format = tmap_get_reads_file_format_string(opt->_name); \
      tmap_file_fprintf(tmap_file_stderr, "[%s]", reads_format); \
      free(reads_format); \
  }

/*
 * Define the print functions for each opt.
 */
// global options
__tmap_map_opt_option_print_func_chars_init(fn_fasta, "not using")
__tmap_map_opt_option_print_func_char_array_init(fn_reads, fn_reads_num, "stdin")
__tmap_map_opt_option_print_func_reads_format_init(reads_format)
__tmap_map_opt_option_print_func_chars_init(fn_sam, "stdout")
__tmap_map_opt_option_print_func_int_init(score_match)
__tmap_map_opt_option_print_func_int_init(pen_mm)
__tmap_map_opt_option_print_func_int_init(pen_gapo)
__tmap_map_opt_option_print_func_int_init(pen_gape)
__tmap_map_opt_option_print_func_int_init(bw)
__tmap_map_opt_option_print_func_int_init(softclip_type)
__tmap_map_opt_option_print_func_int_init(dup_window)
__tmap_map_opt_option_print_func_int_init(max_seed_band)
__tmap_map_opt_option_print_func_tf_init(unroll_banding)
__tmap_map_opt_option_print_func_int_init(score_thr)
__tmap_map_opt_option_print_func_int_init(reads_queue_size)
__tmap_map_opt_option_print_func_int_autodetected_init(num_threads, num_threads_autodetected)
__tmap_map_opt_option_print_func_int_init(aln_output_mode)
__tmap_map_opt_option_print_func_char_array_init(sam_rg, sam_rg_num, "not using")
__tmap_map_opt_option_print_func_tf_init(bidirectional)
__tmap_map_opt_option_print_func_tf_init(seq_eq)
__tmap_map_opt_option_print_func_tf_init(ignore_rg_sam_tags)
__tmap_map_opt_option_print_func_compr_init(input_compr_gz, input_compr, TMAP_FILE_GZ_COMPRESSION)
__tmap_map_opt_option_print_func_compr_init(input_compr_bz2, input_compr, TMAP_FILE_BZ2_COMPRESSION)
__tmap_map_opt_option_print_func_int_init(shm_key)
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
__tmap_map_opt_option_print_func_double_init(sample_reads)
#endif
__tmap_map_opt_option_print_func_int_init(vsw_type)
__tmap_map_opt_option_print_func_verbosity_init()
// flowspace
__tmap_map_opt_option_print_func_int_init(fscore)
__tmap_map_opt_option_print_func_tf_init(softclip_key)
__tmap_map_opt_option_print_func_tf_init(sam_flowspace_tags)
__tmap_map_opt_option_print_func_tf_init(ignore_flowgram)
__tmap_map_opt_option_print_func_tf_init(aln_flowspace)
// pairing
__tmap_map_opt_option_print_func_int_init(pairing)
__tmap_map_opt_option_print_func_int_init(strandedness)
__tmap_map_opt_option_print_func_int_init(positioning)
__tmap_map_opt_option_print_func_double_init(ins_size_mean)
__tmap_map_opt_option_print_func_double_init(ins_size_std)
__tmap_map_opt_option_print_func_double_init(ins_size_std_max_num)
__tmap_map_opt_option_print_func_tf_init(read_rescue)
__tmap_map_opt_option_print_func_double_init(read_rescue_std_num)
__tmap_map_opt_option_print_func_int_init(read_rescue_mapq_thr)
// map1/map2/map3 options, but specific to each
__tmap_map_opt_option_print_func_int_init(min_seq_len)
__tmap_map_opt_option_print_func_int_init(max_seq_len)
// map1/map3 options
__tmap_map_opt_option_print_func_int_init(seed_length)
// map2/map3 options
__tmap_map_opt_option_print_func_int_init(max_seed_hits)
// map3/map4 options
__tmap_map_opt_option_print_func_double_init(hit_frac)
__tmap_map_opt_option_print_func_int_init(seed_step)
// map1 options
__tmap_map_opt_option_print_func_int_init(seed_max_diff)
__tmap_map_opt_option_print_func_int_init(seed2_length)
__tmap_map_opt_option_print_func_np_init(max_diff, max_diff_fnr)
__tmap_map_opt_option_print_func_double_init(max_err_rate)
__tmap_map_opt_option_print_func_nf_init(max_mm, max_mm_frac)
__tmap_map_opt_option_print_func_nf_init(max_gapo, max_gapo_frac)
__tmap_map_opt_option_print_func_nf_init(max_gape, max_gape_frac)
__tmap_map_opt_option_print_func_int_init(max_cals_del)
__tmap_map_opt_option_print_func_int_init(indel_ends_bound)
__tmap_map_opt_option_print_func_int_init(max_best_cals)
__tmap_map_opt_option_print_func_int_init(max_entries)
// map2 options
__tmap_map_opt_option_print_func_double_init(length_coef)
__tmap_map_opt_option_print_func_int_init(max_seed_intv)
__tmap_map_opt_option_print_func_int_init(z_best)
__tmap_map_opt_option_print_func_int_init(seeds_rev)
__tmap_map_opt_option_print_func_tf_init(narrow_rmdup)
// map3 options
__tmap_map_opt_option_print_func_int_init(hp_diff)
__tmap_map_opt_option_print_func_tf_init(fwd_search)
__tmap_map_opt_option_print_func_double_init(skip_seed_frac)
// map4 options
__tmap_map_opt_option_print_func_int_init(min_seed_length)
__tmap_map_opt_option_print_func_int_init(max_seed_length)
__tmap_map_opt_option_print_func_int_init(max_iwidth)
__tmap_map_opt_option_print_func_int_init(max_repr)
__tmap_map_opt_option_print_func_tf_init(rand_repr)
__tmap_map_opt_option_print_func_tf_init(use_min)
// mapvsw options
// stage options
__tmap_map_opt_option_print_func_int_init(stage_score_thr)
__tmap_map_opt_option_print_func_int_init(stage_mapq_thr)
__tmap_map_opt_option_print_func_tf_init(stage_keep_all)
__tmap_map_opt_option_print_func_double_init(stage_seed_freqc)
__tmap_map_opt_option_print_func_double_init(stage_seed_freqc_group_frac)
__tmap_map_opt_option_print_func_int_init(stage_seed_freqc_rand_repr)
__tmap_map_opt_option_print_func_int_init(stage_seed_freqc_min_groups)
__tmap_map_opt_option_print_func_int_init(stage_seed_max_length)

static int32_t
tmap_map_opt_option_flag_length(tmap_map_opt_option_t *opt)
{
  int32_t flag_length = 0;
  if(0 < opt->option.val) flag_length += 3; // dash, char, comma
  if(NULL != opt->name) flag_length += 2 + strlen(opt->name);
  return flag_length;
}

static int32_t 
tmap_map_opt_option_type_length(tmap_map_opt_option_t *opt)
{
  int32_t type_length = 0;
  if(TMAP_MAP_OPT_TYPE_NONE != opt->type) {
      type_length = strlen(tmap_map_opt_input_types[opt->type]);
  }
  return type_length;
}
          
static void
tmap_map_opt_option_print(tmap_map_opt_option_t *opt, tmap_map_opt_t *parent_opt)
{
  int32_t i, j, flag_length, type_length, length_to_description = 0;
  static char *spacer = "         ";

  flag_length = tmap_map_opt_option_flag_length(opt);
  type_length = tmap_map_opt_option_type_length(opt);

  if(NULL == opt->option.name) {
      tmap_error("option did not have a name", Exit, OutOfRange);
  }

  // spacer
  length_to_description += tmap_file_fprintf(tmap_file_stderr, spacer);
  // short flag, if available
  if(0 < opt->option.val) {
      length_to_description += tmap_file_fprintf(tmap_file_stderr, "-%c,", (char)opt->option.val);
  }
  // long flag
  length_to_description += tmap_file_fprintf(tmap_file_stderr, "--%s", opt->option.name);
  if(NULL != parent_opt) {
      for(i=flag_length;i< parent_opt->options->max_flag_length;i++) {
          length_to_description += tmap_file_fprintf(tmap_file_stderr, " ");
      }
  }
  length_to_description += tmap_file_fprintf(tmap_file_stderr, " ");
  // type
  length_to_description += tmap_file_fprintf(tmap_file_stderr, "%s",
                    (TMAP_MAP_OPT_TYPE_NONE == opt->type) ? "" : tmap_map_opt_input_types[opt->type]);
  if(NULL != parent_opt) {
      for(i=type_length;i<parent_opt->options->max_type_length;i++) {
          length_to_description += tmap_file_fprintf(tmap_file_stderr, " ");
      }
  }
  // spacer
  length_to_description += tmap_file_fprintf(tmap_file_stderr, spacer);
  // description
  tmap_file_fprintf(tmap_file_stderr, "%s",
                    opt->description);
  // value, if available
  if(NULL != opt->print_func && NULL != parent_opt) {
      tmap_file_fprintf(tmap_file_stderr, " ");
      opt->print_func(parent_opt);
  }
  tmap_file_fprintf(tmap_file_stderr, "\n");
  // multi-option description
  if(NULL != opt->multi_options) {
      for(i=0;NULL != opt->multi_options[i];i++) {
          // spacers
          for(j=0;j<length_to_description;j++) {
              tmap_file_fprintf(tmap_file_stderr, " ");
          }
          tmap_file_fprintf(tmap_file_stderr, spacer);
          tmap_file_fprintf(tmap_file_stderr, "%s\n", opt->multi_options[i]);
      }
  }
}

static tmap_map_opt_options_t *
tmap_map_opt_options_init()
{
  return tmap_calloc(1, sizeof(tmap_map_opt_options_t), "");
}

/*
  tmap_map_opt_options_add(opt->options, "seed-length", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the k-mer length to seed CALs (-1 to disable)",
                           NULL,
                           tmap_map_opt_option_print_func_seed_length,
                           TMAP_MAP_ALGO_MAP1 | TMAP_MAP_ALGO_MAP3);
*/

static void
tmap_map_opt_options_add(tmap_map_opt_options_t *options, const char *name, 
                         int32_t has_arg, int32_t *flag, int32_t val, int32_t type, 
                         char *description, char **multi_options,
                         tmap_map_opt_option_print_t print_func,
                         int32_t algos)
{
  int32_t flag_length = 0, type_length = 0;
  while(options->mem <= options->n) {
      if(0 == options->mem) {
          options->mem = 4;
      }
      else {
          options->mem *= 2;
      }
      options->options = tmap_realloc(options->options, sizeof(tmap_map_opt_option_t) * options->mem, "options->options");
  }
  
  if(NULL == name) {
      tmap_error("option did not have a name", Exit, OutOfRange);
  }
  if(NULL == description) {
      tmap_error("option did not have a description", Exit, OutOfRange);
  }
  options->options[options->n].name = tmap_strdup(name);
  options->options[options->n].option.name = options->options[options->n].name;
  options->options[options->n].option.has_arg = has_arg;
  options->options[options->n].option.flag = flag;
  options->options[options->n].option.val = val;
  options->options[options->n].type = type;
  options->options[options->n].description = (NULL == description) ? NULL : tmap_strdup(description);
  options->options[options->n].multi_options = multi_options;
  options->options[options->n].print_func = (NULL == print_func) ? NULL : print_func;
  options->options[options->n].algos = algos;

  flag_length = tmap_map_opt_option_flag_length(&options->options[options->n]);
  if(options->max_flag_length < flag_length) {
      options->max_flag_length = flag_length;
  }

  type_length = tmap_map_opt_option_type_length(&options->options[options->n]);
  if(options->max_type_length < type_length) {
      options->max_type_length = type_length;
  }

  options->n++;
}

static void
tmap_map_opt_options_destroy(tmap_map_opt_options_t *options)
{
  int32_t i;
  for(i=0;i<options->n;i++) {
      free(options->options[i].name);
      free(options->options[i].description);
  }
  free(options->options);
  free(options);
}

static void
tmap_map_opt_init_helper(tmap_map_opt_t *opt)
{
  static char *softclipping_type[] = {"0 - allow on the left and right portion of the read",
      "1 - allow on the left portion of the read",
      "2 - allow on the right portion of the read",
      "3 - do not allow soft-clipping",
      NULL};
  static char *aln_output_mode[] = {"0 - unique best hits",
      "1 - random best hit",
      "2 - all best hits",
      "3 - all alignments",
      NULL};
  static char *vsw_type[] = {
      "NB: currently only #1, #4, and #6 have been tested",
      "1 - lh3/ksw.c/nh13",
      "2 - simple VSW",
      "3 - SHRiMP2 VSW [not working]",
      "4 - Psyho (Top Coder #1)",
      "5 - ACRush (Top Coder #2)",
      "6 - folsena (Top Coder #3)",
      "7 - logicmachine (Top Coder #4)",
      "8 - venco (Top Coder #5) [not working]",
      "9 - Bladze (Top Coder #6)",
      "10 - ngthuydiem (Top Coder #7) [Farrar cut-and-paste]",
      NULL};
  static char *pairing[] = {"0 - no pairing is to be performed", "1 - mate pairs (-S 0 -P 1)", "2 - paired end (-S 1 -P 0)", NULL};
  static char *strandedness[] = {"0 - same strand", "1 - opposite strand", NULL};
  static char *positioning[] = {"0 - read one before read two", "1 - read two before read one", NULL};

  opt->options = tmap_map_opt_options_init();

  // global options
  tmap_map_opt_options_add(opt->options, "fn-fasta", required_argument, 0, 'f', 
                           TMAP_MAP_OPT_TYPE_FILE,
                           "FASTA reference file name",
                           NULL,
                           tmap_map_opt_option_print_func_fn_fasta,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "fn-reads", required_argument, 0, 'r', 
                           TMAP_MAP_OPT_TYPE_FILE,
                           "the reads file name", 
                           NULL,
                           tmap_map_opt_option_print_func_fn_reads,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "reads-format", required_argument, 0, 'i', 
                           TMAP_MAP_OPT_TYPE_STRING,
                           "the reads file format (opt->options, fastq|fq|fasta|fa|sff|sam|bam)",
                           NULL,
                           tmap_map_opt_option_print_func_reads_format,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "fn-sam", required_argument, 0, 's', 
                           TMAP_MAP_OPT_TYPE_FILE,
                           "the SAM file name",
                           NULL,
                           tmap_map_opt_option_print_func_fn_sam,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "score-match", required_argument, 0, 'A', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "score for a match",
                           NULL,
                           tmap_map_opt_option_print_func_score_match,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "pen-mismatch", required_argument, 0, 'M', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the mismatch penalty",
                           NULL,
                           tmap_map_opt_option_print_func_pen_mm,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "pen-gap-open", required_argument, 0, 'O', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the indel start penalty",
                           NULL,
                           tmap_map_opt_option_print_func_pen_gapo,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "pen-gap-extension", required_argument, 0, 'E', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the indel extension penalty",
                           NULL,
                           tmap_map_opt_option_print_func_pen_gape,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "band-width", required_argument, 0, 'w', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the band width",
                           NULL,
                           tmap_map_opt_option_print_func_bw,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "softclip-type", required_argument, 0, 'g', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the soft-clipping type",
                           softclipping_type,
                           tmap_map_opt_option_print_func_softclip_type,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "duplicate-window", required_argument, 0, 'W', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "remove duplicate alignments within this bp window (-1 to disable)",
                           NULL,
                           tmap_map_opt_option_print_func_dup_window,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "max-seed-band", required_argument, 0, 'B', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the window of bases in which to group seeds",
                           NULL,
                           tmap_map_opt_option_print_func_max_seed_band,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "unroll-banding", no_argument, 0, 'U', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "unroll the grouped seeds from banding if multiple alignments are found",
                           NULL,
                           tmap_map_opt_option_print_func_unroll_banding,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "score-thres", required_argument, 0, 'T', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "score threshold divided by the match score",
                           NULL,
                           tmap_map_opt_option_print_func_score_thr,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "reads-queue-size", required_argument, 0, 'q', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the queue size for the reads (-1 to disable)",
                           NULL,
                           tmap_map_opt_option_print_func_reads_queue_size,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "num-threads", required_argument, 0, 'n', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the number of threads",
                           NULL,
                           tmap_map_opt_option_print_func_num_threads,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "aln-output-mode", required_argument, 0, 'a', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "output filter",
                           aln_output_mode,
                           tmap_map_opt_option_print_func_aln_output_mode,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "sam-read-group", required_argument, 0, 'R', 
                           TMAP_MAP_OPT_TYPE_STRING,
                           "the RG tags to add to the SAM header",
                           NULL,
                           tmap_map_opt_option_print_func_sam_rg,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "bidirectional", no_argument, 0, 'D', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "specifies the input reads are to be annotated as bidirectional",
                           NULL,
                           tmap_map_opt_option_print_func_bidirectional,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "use-seq-equal", no_argument, 0, 'I', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "specifies to use the '=' symbol in the SEQ field",
                           NULL,
                           tmap_map_opt_option_print_func_seq_eq,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "ignore-rg-from-sam", no_argument, 0, 'C', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "specifies to not use the RG header and RG record tags in the SAM file",
                           NULL,
                           tmap_map_opt_option_print_func_ignore_rg_sam_tags,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "input-gz", no_argument, 0, 'z', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "the input is gz (gzip) compressed",
                           NULL,
                           tmap_map_opt_option_print_func_input_compr_gz,
                           TMAP_MAP_ALGO_GLOBAL);
#ifndef DISABLE_BZ2
  tmap_map_opt_options_add(opt->options, "input-bz2", no_argument, 0, 'j', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "the input is bz2 (bzip2) compressed",
                           NULL,
                           tmap_map_opt_option_print_func_input_compr_bz2,
                           TMAP_MAP_ALGO_GLOBAL);
#endif
  tmap_map_opt_options_add(opt->options, "shared-memory-key", required_argument, 0, 'k', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "use shared memory with the following key",
                           NULL,
                           tmap_map_opt_option_print_func_shm_key,
                           TMAP_MAP_ALGO_GLOBAL);
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
  tmap_map_opt_options_add(opt->options, "sample-reads", required_argument, 0, 'x',
                           TMAP_MAP_OPT_TYPE_FLOAT,
                           "sample the reads at this fraction",
                           NULL,
                           tmap_map_opt_option_print_func_sample_reads,
                           TMAP_MAP_ALGO_GLOBAL);
#endif
  tmap_map_opt_options_add(opt->options, "vsw-type", required_argument, 0, 'H',
                           TMAP_MAP_OPT_TYPE_INT,
                           "the vectorized smith-waterman algorithm (very untested)",
                           vsw_type,
                           tmap_map_opt_option_print_func_vsw_type,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "help", no_argument, 0, 'h', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "print this message",
                           NULL,
                           NULL,
                           TMAP_MAP_ALGO_GLOBAL);
  tmap_map_opt_options_add(opt->options, "verbose", no_argument, 0, 'v', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "print verbose progress information",
                           NULL,
                           tmap_map_opt_option_print_func_verbosity,
                           TMAP_MAP_ALGO_GLOBAL);
  
  // flowspace options
  tmap_map_opt_options_add(opt->options, "pen-flow-error", required_argument, 0, 'X', 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the flow score penalty",
                           NULL,
                           tmap_map_opt_option_print_func_fscore,
                           TMAP_MAP_ALGO_FLOWSPACE);
  tmap_map_opt_options_add(opt->options, "softclip-key", no_argument, 0, 'y', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "soft clip only the last base of the key",
                           NULL,
                           tmap_map_opt_option_print_func_softclip_key,
                           TMAP_MAP_ALGO_FLOWSPACE);
  tmap_map_opt_options_add(opt->options, "sam-flowspace-tags", no_argument, 0, 'Y', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "include flow space specific SAM tags when available",
                           NULL,
                           tmap_map_opt_option_print_func_sam_flowspace_tags,
                           TMAP_MAP_ALGO_FLOWSPACE);
  tmap_map_opt_options_add(opt->options, "ignore-flowgram", no_argument, 0, 'N', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "do not use the flowgram, otherwise use the flowgram when available",
                           NULL,
                           tmap_map_opt_option_print_func_ignore_flowgram,
                           TMAP_MAP_ALGO_FLOWSPACE);
  tmap_map_opt_options_add(opt->options, "final-flowspace", no_argument, 0, 'F', 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "produce the final alignment in flow space",
                           NULL,
                           tmap_map_opt_option_print_func_aln_flowspace,
                           TMAP_MAP_ALGO_FLOWSPACE);

  // pairing options
  tmap_map_opt_options_add(opt->options, "pairing", required_argument, 0, 'Q',
                           TMAP_MAP_OPT_TYPE_INT,
                           "the insert pairing",
                           pairing,
                           tmap_map_opt_option_print_func_pairing,
                           TMAP_MAP_ALGO_PAIRING);
  tmap_map_opt_options_add(opt->options, "strandedness", required_argument, 0, 'S',
                           TMAP_MAP_OPT_TYPE_INT,
                           "the insert strandedness",
                           strandedness,
                           tmap_map_opt_option_print_func_strandedness,
                           TMAP_MAP_ALGO_PAIRING);
  tmap_map_opt_options_add(opt->options, "positioning", required_argument, 0, 'P',
                           TMAP_MAP_OPT_TYPE_INT,
                           "the insert positioning when on the same strand (-S 0)",
                           positioning,
                           tmap_map_opt_option_print_func_positioning,
                           TMAP_MAP_ALGO_PAIRING);
  tmap_map_opt_options_add(opt->options, "ins-size-mean", required_argument, 0, 'b',
                           TMAP_MAP_OPT_TYPE_FLOAT,
                           "the mean insert size",
                           NULL,
                           tmap_map_opt_option_print_func_ins_size_mean,
                           TMAP_MAP_ALGO_PAIRING);
  tmap_map_opt_options_add(opt->options, "ins-size-std", required_argument, 0, 'c',
                           TMAP_MAP_OPT_TYPE_FLOAT,
                           "the insert size standard deviation",
                           NULL,
                           tmap_map_opt_option_print_func_ins_size_std,
                           TMAP_MAP_ALGO_PAIRING);
  tmap_map_opt_options_add(opt->options, "ins-size-std-max-num", required_argument, 0, 'd',
                           TMAP_MAP_OPT_TYPE_FLOAT,
                           "the insert size maximum standard deviation",
                           NULL,
                           tmap_map_opt_option_print_func_ins_size_std_max_num,
                           TMAP_MAP_ALGO_PAIRING);
  tmap_map_opt_options_add(opt->options, "read-rescue", no_argument, 0, 'L',
                           TMAP_MAP_OPT_TYPE_NONE,
                           "perform read rescue",
                           NULL,
                           tmap_map_opt_option_print_func_read_rescue,
                           TMAP_MAP_ALGO_PAIRING);
  tmap_map_opt_options_add(opt->options, "read-rescue-std-num", required_argument, 0, 'l',
                           TMAP_MAP_OPT_TYPE_FLOAT,
                           "the number of standard deviations around the mean insert size to perform read rescue",
                           NULL,
                           tmap_map_opt_option_print_func_read_rescue_std_num,
                           TMAP_MAP_ALGO_PAIRING);
  tmap_map_opt_options_add(opt->options, "read-rescue-mapq-thr", required_argument, 0, 'm',
                           TMAP_MAP_OPT_TYPE_INT,
                           "mapping quality threshold for read rescue",
                           NULL,
                           tmap_map_opt_option_print_func_read_rescue_mapq_thr,
                           TMAP_MAP_ALGO_PAIRING);

  // map1/map3 options
  tmap_map_opt_options_add(opt->options, "seed-length", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the k-mer length to seed CALs (-1 to disable)",
                           NULL,
                           tmap_map_opt_option_print_func_seed_length,
                           TMAP_MAP_ALGO_MAP1 | TMAP_MAP_ALGO_MAP3);

  // map2/map3 options
  tmap_map_opt_options_add(opt->options, "max-seed-hits", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the maximum number of hits returned by a seed",
                           NULL,
                           tmap_map_opt_option_print_func_max_seed_hits,
                           TMAP_MAP_ALGO_MAP2 | TMAP_MAP_ALGO_MAP3);

  // map3/map4 options
  tmap_map_opt_options_add(opt->options, "hit-frac", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_FLOAT,
                           "the fraction of seed positions that are under the maximum",
                           NULL,
                           tmap_map_opt_option_print_func_hit_frac,
                           TMAP_MAP_ALGO_MAP3 | TMAP_MAP_ALGO_MAP4);
  tmap_map_opt_options_add(opt->options, "seed-step", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the number of bases to increase the seed for each seed increase iteration (-1 to disable)",
                           NULL,
                           tmap_map_opt_option_print_func_seed_step,
                           TMAP_MAP_ALGO_MAP3 | TMAP_MAP_ALGO_MAP4);

  // map1 options
  tmap_map_opt_options_add(opt->options, "seed-max-diff", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "maximum number of edits in the seed",
                           NULL,
                           tmap_map_opt_option_print_func_seed_max_diff,
                           TMAP_MAP_ALGO_MAP1);
  tmap_map_opt_options_add(opt->options, "seed2-length", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the secondary seed length (-1 to disable)",
                           NULL,
                           tmap_map_opt_option_print_func_seed2_length,
                           TMAP_MAP_ALGO_MAP1);
  tmap_map_opt_options_add(opt->options, "max-diff", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_NUM,
                           "maximum number of edits or false-negative probability assuming the maximum error rate",
                           NULL,
                           tmap_map_opt_option_print_func_max_diff,
                           TMAP_MAP_ALGO_MAP1);
  tmap_map_opt_options_add(opt->options, "max-error-rate", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_FLOAT,
                           "the assumed per-base maximum error rate",
                           NULL,
                           tmap_map_opt_option_print_func_max_err_rate,
                           TMAP_MAP_ALGO_MAP1);
  tmap_map_opt_options_add(opt->options, "max-mismatches", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_NUM,
                           "maximum number of or (read length) fraction of mismatches",
                           NULL,
                           tmap_map_opt_option_print_func_max_mm,
                           TMAP_MAP_ALGO_MAP1);
  tmap_map_opt_options_add(opt->options, "max-gap-opens", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_NUM,
                           "maximum number of or (read length) fraction of indel starts",
                           NULL,
                           tmap_map_opt_option_print_func_max_gapo,
                           TMAP_MAP_ALGO_MAP1);
  tmap_map_opt_options_add(opt->options, "max-gap-extensions", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_NUM,
                           "maximum number of or (read length) fraction of indel extensions",
                           NULL,
                           tmap_map_opt_option_print_func_max_gape,
                           TMAP_MAP_ALGO_MAP1);
  tmap_map_opt_options_add(opt->options, "max-cals-deletion", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the maximum number of CALs to extend a deletion ",
                           NULL,
                           tmap_map_opt_option_print_func_max_cals_del,
                           TMAP_MAP_ALGO_MAP1);
  tmap_map_opt_options_add(opt->options, "indel-ends-bound", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "number of bps from the end of the read ",
                           NULL,
                           tmap_map_opt_option_print_func_indel_ends_bound,
                           TMAP_MAP_ALGO_MAP1);
  tmap_map_opt_options_add(opt->options, "max-best-cals", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "optimal CALs have been found ",
                           NULL,
                           tmap_map_opt_option_print_func_max_best_cals,
                           TMAP_MAP_ALGO_MAP1);
  tmap_map_opt_options_add(opt->options, "max-nodes", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "maximum number of alignment nodes ",
                           NULL,
                           tmap_map_opt_option_print_func_max_entries,
                           TMAP_MAP_ALGO_MAP1);

  // map2 options
  tmap_map_opt_options_add(opt->options, "length-coef", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_FLOAT,
"coefficient of length-threshold adjustment",
                           NULL,
                           tmap_map_opt_option_print_func_length_coef,
                           TMAP_MAP_ALGO_MAP2);
  tmap_map_opt_options_add(opt->options, "max-seed-intv", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
"maximum seeding interval size",
                           NULL,
                           tmap_map_opt_option_print_func_max_seed_intv,
                           TMAP_MAP_ALGO_MAP2);
  tmap_map_opt_options_add(opt->options, "z-best", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the maximum number of top-scoring nodes to keep on each iteration",
                           NULL,
                           tmap_map_opt_option_print_func_z_best,
                           TMAP_MAP_ALGO_MAP2);
  tmap_map_opt_options_add(opt->options, "seeds-rev", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
"# seeds to trigger reverse alignment",
                           NULL,
                           tmap_map_opt_option_print_func_seeds_rev,
                           TMAP_MAP_ALGO_MAP2);
  tmap_map_opt_options_add(opt->options, "narrow-rmdup", no_argument, 0, 0,
                           TMAP_MAP_OPT_TYPE_NONE,
                           "remove duplicates for narrow SA hits",
                           NULL,
                           tmap_map_opt_option_print_func_narrow_rmdup,
                           TMAP_MAP_ALGO_MAP2);

  // map3 options
  tmap_map_opt_options_add(opt->options, "hp-diff", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "single homopolymer error difference for enumeration",
                           NULL,
                           tmap_map_opt_option_print_func_hp_diff,
                           TMAP_MAP_ALGO_MAP3);
  tmap_map_opt_options_add(opt->options, "fwd-search", no_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "use forward search instead of a reverse search",
                           NULL,
                           tmap_map_opt_option_print_func_fwd_search,
                           TMAP_MAP_ALGO_MAP3);
  tmap_map_opt_options_add(opt->options, "skip-seed-frac", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_FLOAT,
                           "the fraction of a seed to skip when a lookup succeeds",
                           NULL,
                           tmap_map_opt_option_print_func_skip_seed_frac,
                           TMAP_MAP_ALGO_MAP3);
  // map4
  tmap_map_opt_options_add(opt->options, "min-seed-length", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the minimum seed length to accept hits",
                           NULL,
                           tmap_map_opt_option_print_func_min_seed_length,
                           TMAP_MAP_ALGO_MAP4);
  tmap_map_opt_options_add(opt->options, "max-seed-length", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the maximum seed length to accept hits",
                           NULL,
                           tmap_map_opt_option_print_func_max_seed_length,
                           TMAP_MAP_ALGO_MAP4);
  tmap_map_opt_options_add(opt->options, "max-iwidth", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the maximum interval size to accept a hit",
                           NULL,
                           tmap_map_opt_option_print_func_max_iwidth,
                           TMAP_MAP_ALGO_MAP4);
  tmap_map_opt_options_add(opt->options, "max-repr", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the maximum representitive hits for repetitive hits",
                           NULL,
                           tmap_map_opt_option_print_func_max_repr,
                           TMAP_MAP_ALGO_MAP4);
  tmap_map_opt_options_add(opt->options, "rand-repr", no_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "choose the representitive hits randomly, otherwise uniformly",
                           NULL,
                           tmap_map_opt_option_print_func_rand_repr,
                           TMAP_MAP_ALGO_MAP4);
  tmap_map_opt_options_add(opt->options, "use-min", no_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "when seed stepping, try seeding when at least the minimum seed length is present, otherwise maximum",
                           NULL,
                           tmap_map_opt_option_print_func_use_min,
                           TMAP_MAP_ALGO_MAP4);

  // mapvsw options
  // None

  // map1/map2/map3 options, but specific to each
  tmap_map_opt_options_add(opt->options, "min-seq-length", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the minimum sequence length to examine (-1 to disable)",
                           NULL,
                           tmap_map_opt_option_print_func_min_seq_len,
                           ~(TMAP_MAP_ALGO_MAPALL | TMAP_MAP_ALGO_STAGE));
  tmap_map_opt_options_add(opt->options, "max-seq-length", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the maximum sequence length to examine",
                           NULL,
                           tmap_map_opt_option_print_func_max_seq_len,
                           ~(TMAP_MAP_ALGO_MAPALL | TMAP_MAP_ALGO_STAGE));

  // stage options
  tmap_map_opt_options_add(opt->options, "stage-score-thres", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "score threshold for the stage divided by the match score",
                           NULL,
                           tmap_map_opt_option_print_func_stage_score_thr,
                           TMAP_MAP_ALGO_STAGE);
  tmap_map_opt_options_add(opt->options, "stage-mapq-thres", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "mapping quality threshold for the stage divided by the match score",
                           NULL,
                           tmap_map_opt_option_print_func_stage_mapq_thr,
                           TMAP_MAP_ALGO_STAGE);
  tmap_map_opt_options_add(opt->options, "stage-keep-all", no_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_NONE,
                           "keep mappings from the first stage for the next stage",
                           NULL,
                           tmap_map_opt_option_print_func_stage_keep_all,
                           TMAP_MAP_ALGO_STAGE);
  tmap_map_opt_options_add(opt->options, "stage-seed-freq-cutoff", required_argument, 0, 0,
                           TMAP_MAP_OPT_TYPE_FLOAT,
                           "the minimum frequency of a seed to be considered for mapping",
                           NULL,
                           tmap_map_opt_option_print_func_stage_seed_freqc,
                           TMAP_MAP_ALGO_STAGE);
  tmap_map_opt_options_add(opt->options, "stage-seed-freq-cutoff-group-frac", required_argument, 0, 0,
                           TMAP_MAP_OPT_TYPE_FLOAT,
                           "if more than this fraction of groups were filtered, keep representative hits",
                           NULL,
                           tmap_map_opt_option_print_func_stage_seed_freqc_group_frac,
                           TMAP_MAP_ALGO_STAGE);
  tmap_map_opt_options_add(opt->options, "stage-seed-freq-cutoff-rand-repr", required_argument, 0, 0,
                           TMAP_MAP_OPT_TYPE_INT,
                           "the number of representative hits to keep",
                           NULL,
                           tmap_map_opt_option_print_func_stage_seed_freqc_rand_repr,
                           TMAP_MAP_ALGO_STAGE);
  tmap_map_opt_options_add(opt->options, "stage-seed-freq-cutoff-min-groups", required_argument, 0, 0,
                           TMAP_MAP_OPT_TYPE_INT,
                           "the minimum of groups required after the filter has been applied, otherwise iteratively reduce the filter",
                           NULL,
                           tmap_map_opt_option_print_func_stage_seed_freqc_min_groups,
                           TMAP_MAP_ALGO_STAGE);
  tmap_map_opt_options_add(opt->options, "stage-seed-max-length", required_argument, 0, 0, 
                           TMAP_MAP_OPT_TYPE_INT,
                           "the length of the prefix of the read to consider during seeding",
                           NULL,
                           tmap_map_opt_option_print_func_stage_seed_max_length,
                           TMAP_MAP_ALGO_STAGE);

  /*
  // Prints out all single-flag command line options
  int i, c;
  for(c=1;c<256;c++) {
      for(i=0;i<opt->options->n;i++) {
          if(c == opt->options->options[i].option.val) {
              fputc((char)c, stderr);
              if(required_argument == opt->options->options[i].option.has_arg) {
                  fputc(':', stderr);
              }
          }
      }
  }
  fprintf(stderr, "\n");
  */
  
  /*
  // Prints out all command line options
  int i, c;
  for(c=0;c<256;c++) {
      for(i=0;i<opt->options->n;i++) {
          if(c == opt->options->options[i].option.val) {
              if(0 < opt->options->options[i].option.val) {
                  tmap_file_fprintf(tmap_file_stderr, "-%c,", (char)opt->options->options[i].option.val);
              }
              // long flag
              tmap_file_fprintf(tmap_file_stderr, "--%s\n", opt->options->options[i].option.name);
          }
      }
  }
  fprintf(stderr, "\n");
  */
}

tmap_map_opt_t *
tmap_map_opt_init(int32_t algo_id)
{
  tmap_map_opt_t *opt = NULL;

  opt = tmap_calloc(1, sizeof(tmap_map_opt_t), "opt");

  // internal data
  opt->algo_id = algo_id;
  opt->algo_stage = -1;
  opt->argv = NULL;
  opt->argc = -1;

  // global options 
  opt->fn_fasta = NULL;
  opt->fn_reads = NULL;
  opt->fn_reads_num = 0;
  opt->reads_format = TMAP_READS_FORMAT_UNKNOWN;
  opt->fn_sam = NULL;
  opt->score_match = TMAP_MAP_OPT_SCORE_MATCH;
  opt->pen_mm = TMAP_MAP_OPT_PEN_MM;
  opt->pen_gapo = TMAP_MAP_OPT_PEN_GAPO;
  opt->pen_gape = TMAP_MAP_OPT_PEN_GAPE;
  opt->bw = 50; 
  opt->softclip_type = TMAP_MAP_OPT_SOFT_CLIP_RIGHT;
  opt->dup_window = 128;
  opt->max_seed_band = 15;
  opt->unroll_banding = 0;
  opt->score_thr = 8;
  opt->reads_queue_size = 262144;
  opt->num_threads = tmap_detect_cpus();
  opt->num_threads_autodetected = 1;
  opt->aln_output_mode = TMAP_MAP_OPT_ALN_MODE_RAND_BEST;
  opt->sam_rg = NULL;
  opt->sam_rg_num = 0;
  opt->bidirectional = 0;
  opt->seq_eq = 0;
  opt->ignore_rg_sam_tags = 0;
  opt->input_compr = TMAP_FILE_NO_COMPRESSION;
  opt->shm_key = 0;
  opt->min_seq_len = -1;
  opt->max_seq_len = -1;
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
  opt->sample_reads = 1.0;
#endif
  opt->vsw_type = 4;

  // flowspace options
  opt->fscore = TMAP_MAP_OPT_FSCORE;
  opt->softclip_key = 0;
  opt->sam_flowspace_tags = 0;
  opt->ignore_flowgram = 0;
  opt->aln_flowspace = 0;

  // pairing options
  opt->pairing = 0;
  opt->strandedness = -1;
  opt->positioning = -1;
  opt->ins_size_mean = -1.0;
  opt->ins_size_std = -1.0;
  opt->ins_size_std_max_num  = -1.0;
  opt->read_rescue = 0;
  opt->read_rescue_std_num = -1.0;
  opt->read_rescue_mapq_thr = 0;

  switch(algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      // map1
      opt->seed_length = 32;
      opt->seed_length_set = 1;
      opt->seed_max_diff = 2;
      opt->seed2_length = 48;
      opt->max_diff = -1; opt->max_diff_fnr = 0.04;
      opt->max_mm = 3; opt->max_mm_frac = -1.;
      opt->max_err_rate = 0.02;
      opt->max_gapo = 1; opt->max_gapo_frac = -1.;
      opt->max_gape = 6; opt->max_gape_frac = -1.;
      opt->max_cals_del = 10;
      opt->indel_ends_bound = 5;
      opt->max_best_cals = 32;
      opt->max_entries = 2000000;
      break;
    case TMAP_MAP_ALGO_MAP2:
      // map2
      //opt->mask_level = 0.50; 
      opt->length_coef = 5.5f;
      opt->max_seed_intv = 6; 
      opt->z_best = 1; 
      opt->seeds_rev = 5;
      opt->narrow_rmdup = 0;
      opt->max_seed_hits = 1024;
      break;
    case TMAP_MAP_ALGO_MAP3:
      // map3
      opt->seed_length = -1;
      opt->seed_length_set = 0;
      opt->max_seed_hits = 20;
      opt->hp_diff = 0;
      opt->hit_frac = 0.20;
      opt->seed_step = 8;
      opt->fwd_search = 0;
      opt->skip_seed_frac = 0.2;
      break;
    case TMAP_MAP_ALGO_MAP4:
      // map4
      opt->min_seed_length = -1;
      opt->max_seed_length = 48; 
      opt->hit_frac = 0.20;
      opt->seed_step = 8;
      opt->max_iwidth = 20;
      opt->max_repr = 3;
      opt->rand_repr = 0;
      opt->use_min = 0;
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      // mapvsw
      break;
    case TMAP_MAP_ALGO_STAGE:
      // stage
      opt->stage_score_thr = 8;
      opt->stage_mapq_thr = 23; // 0.5% error
      opt->stage_keep_all = 0;
      opt->stage_seed_freqc = 0.0; //all-pass filter as default
      opt->stage_seed_freqc_group_frac = 0.9; 
      opt->stage_seed_freqc_rand_repr = 2; 
      opt->stage_seed_freqc_min_groups = 1; 
      opt->stage_seed_max_length = -1;
      break;
    default:
      break;
  }

  // build options for parsing and printing
  tmap_map_opt_init_helper(opt);

  opt->sub_opts = NULL;
  opt->num_sub_opts = 0;

  return opt;
}

tmap_map_opt_t*
tmap_map_opt_add_sub_opt(tmap_map_opt_t *opt, int32_t algo_id)
{
  opt->num_sub_opts++;
  opt->sub_opts = tmap_realloc(opt->sub_opts, opt->num_sub_opts * sizeof(tmap_map_opt_t*), "opt->sub_opts");
  opt->sub_opts[opt->num_sub_opts-1] = tmap_map_opt_init(algo_id);
  // copy global options
  tmap_map_opt_copy_global(opt->sub_opts[opt->num_sub_opts-1], opt);
  return opt->sub_opts[opt->num_sub_opts-1];
}

void
tmap_map_opt_destroy(tmap_map_opt_t *opt)
{
  int32_t i;

  free(opt->fn_fasta);
  for(i=0;i<opt->fn_reads_num;i++) {
      free(opt->fn_reads[i]); 
  }
  free(opt->fn_reads);
  free(opt->fn_sam);
  for(i=0;i<opt->sam_rg_num;i++) {
      free(opt->sam_rg[i]);
  }
  free(opt->sam_rg);

  for(i=0;i<opt->num_sub_opts;i++) {
      tmap_map_opt_destroy(opt->sub_opts[i]);
  }
  free(opt->sub_opts);

  // destroy options for parsing and printing
  tmap_map_opt_options_destroy(opt->options);

  free(opt);
}

static void
tmap_map_opt_usage_algo(tmap_map_opt_t *opt, int32_t stage)
{
  int32_t i;
  if(opt->algo_id & TMAP_MAP_ALGO_MAPALL) {
      return; // NB: there are no MAPALL specific options
  }
  else if(opt->algo_id & TMAP_MAP_ALGO_STAGE) {
      tmap_file_fprintf(tmap_file_stderr, "\nstage%d options: [stage options] [algorithm [algorithm options]]+\n", stage);
  }
  else if(stage < 0) {
      tmap_file_fprintf(tmap_file_stderr, "\n%s options (optional):\n", tmap_algo_id_to_name(opt->algo_id));
  }
  else {
      tmap_file_fprintf(tmap_file_stderr, "\n%s stage%d options (optional):\n", tmap_algo_id_to_name(opt->algo_id), stage);
  }
  for(i=0;i<opt->options->n;i++) {
      tmap_map_opt_option_t *o = &opt->options->options[i];

      if(0 < (o->algos & opt->algo_id)) {
          tmap_map_opt_option_print(o, opt);
      }
  }
}

int
tmap_map_opt_usage(tmap_map_opt_t *opt)
{
  int32_t i, j, prev_stage;
  
  // print global options
  tmap_file_fprintf(tmap_file_stderr, "\n");
  if(opt->algo_id == TMAP_MAP_ALGO_MAPALL) {
      tmap_file_fprintf(tmap_file_stderr, "\n%s [global options] [flowspace options] [stage[0-9]+ [stage options] [algorithm [algorithm options]]+]+\n", 
                        tmap_algo_id_to_name(opt->algo_id));
  }
  else {
      tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [global options] [flowspace options] [%s options]\n", 
                        PACKAGE, 
                        tmap_algo_id_to_name(opt->algo_id),
                        tmap_algo_id_to_name(opt->algo_id));
  }
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "global options:\n");
  for(i=0;i<opt->options->n;i++) {
      tmap_map_opt_option_t *o = &opt->options->options[i];

      if(o->algos == TMAP_MAP_ALGO_GLOBAL) {
          tmap_map_opt_option_print(o, opt);
      }
  }

  // print flowspace options
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "flowspace options:\n");
  for(i=0;i<opt->options->n;i++) {
      tmap_map_opt_option_t *o = &opt->options->options[i];

      if(o->algos == TMAP_MAP_ALGO_FLOWSPACE) {
          tmap_map_opt_option_print(o, opt);
      }
  }

  // print pairing options
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "pairing options:\n");
  for(i=0;i<opt->options->n;i++) {
      tmap_map_opt_option_t *o = &opt->options->options[i];

      if(o->algos == TMAP_MAP_ALGO_PAIRING) {
          tmap_map_opt_option_print(o, opt);
      }
  }

  // print algorithm specific options
  for(i=0,prev_stage=-1;i<opt->num_sub_opts;i++) {
      // print the stage
      if(opt->sub_opts[i]->algo_stage != prev_stage) {
          prev_stage = opt->sub_opts[i]->algo_stage;
          // print the stage
          //tmap_file_fprintf(tmap_file_stderr, "\nstage%d options:\n", prev_stage);
          for(j=0;j<opt->options->n;j++) {
              tmap_map_opt_option_t *o = &opt->options->options[j];
              if(o->algos == TMAP_MAP_ALGO_STAGE) {
                  tmap_map_opt_option_print(o, opt);
              }
          }
      }
      
      tmap_map_opt_usage_algo(opt->sub_opts[i], opt->sub_opts[i]->algo_stage);
  }
  tmap_map_opt_usage_algo(opt, -1);

  tmap_map_opt_destroy(opt);

  return 1;
}

static void
tmap_map_opt_add_tabbed(char ***dest, int32_t *dest_num, char *str)
{
  int32_t i, j, l;

  l = strlen(str);
  i = j = 0;
  while(i < l) {
      j = i;
      while(j < l && '\t' != str[j]) {
          j++;
      }
      // add
      (*dest_num)++;
      (*dest) = tmap_realloc((*dest), sizeof(char*) * (*dest_num), "(*dest)");
      (*dest)[(*dest_num)-1] = tmap_malloc(sizeof(char) * (j - i + 1), "(*dest)[(*dest_num)-1]");
      strncpy((*dest)[(*dest_num)-1], str + i, (j - i));
      (*dest)[(*dest_num)-1][(j-i)] = '\0';
      tmap_chomp((*dest)[(*dest_num)-1]); // remove trailing white spaces
      i = j + 1;
  }
}

int32_t
tmap_map_opt_parse(int argc, char *argv[], tmap_map_opt_t *opt)
{
  int i, c, option_index, val = 0;
  char *getopt_format = NULL; 
  int32_t getopt_format_mem = 0;
  int32_t getopt_format_len = 1;
  struct option *options = NULL;

  // set options passed in
  opt->argc = argc; opt->argv = argv;
  
  if(argc == optind) {
      // no need to parse
      //fprintf(stderr, "\n[opt_parse] argc==optind.  no need to parse\n");
      return 1;
  }

  /*
  fprintf(stderr, "\nargc: %d optind: %d\n", argc, optind);
  for(i=optind;i<argc;i++) {
      fprintf(stderr, "[opt_parse] i=%d argv[i]=%s\n", i, argv[i]);
  }
  */
  

  // allocate
  options = tmap_calloc(1, sizeof(struct option) * opt->options->n, "options");

  // format
  getopt_format_len = 0;
  getopt_format_mem = 4;
  getopt_format = tmap_calloc(getopt_format_mem, sizeof(char) * getopt_format_mem, "getopt_format"); 

  // shallow copy
  for(i=0;i<opt->options->n;i++) {
      options[i] = opt->options->options[i].option; 
      if(0 != options[i].val) {
          while(getopt_format_mem < getopt_format_len + 4) {
              getopt_format_mem <<= 1;
              getopt_format = tmap_realloc(getopt_format, sizeof(char) * getopt_format_mem, "getopt_format"); 
          }
          getopt_format[getopt_format_len] = (char)(options[i].val);
          getopt_format_len++;
          getopt_format[getopt_format_len] = '\0';
          if(no_argument != options[i].has_arg) {
              getopt_format[getopt_format_len] = ':';
              getopt_format_len++;
          }
      }
  }
  getopt_format[getopt_format_len] = '\0';

  // check algorithm
  switch(opt->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
    case TMAP_MAP_ALGO_MAP2:
    case TMAP_MAP_ALGO_MAP3:
    case TMAP_MAP_ALGO_MAP4:
    case TMAP_MAP_ALGO_MAPVSW:
    case TMAP_MAP_ALGO_STAGE:
    case TMAP_MAP_ALGO_MAPALL:
      break;
    default:
      tmap_error("unrecognized algorithm", Exit, OutOfRange);
      break;
  }
  
  while((c = getopt_long(argc, argv, getopt_format, options, &option_index)) >= 0) {
      // Global options
      if(c == '?') {
          break;
      }
      else if(c == 'A' || (0 == c && 0 == strcmp("score-match", options[option_index].name))) {       
          opt->score_match = atoi(optarg);
      }
      else if(c == 'B' || (0 == c && 0 == strcmp("max-seed-band", options[option_index].name))) {       
          opt->max_seed_band = atoi(optarg);
      }
      else if(c == 'C' || (0 == c && 0 == strcmp("ignore-rg-from-sam", options[option_index].name))) {
          opt->ignore_rg_sam_tags = 1;
      }
      else if(c == 'D' || (0 == c && 0 == strcmp("bidirectional", options[option_index].name))) {       
          opt->bidirectional = 1;
      }
      else if(c == 'U' || (0 == c && 0 == strcmp("unroll-banding", options[option_index].name))) {       
          opt->unroll_banding = 1;
      }
      else if(c == 'E' || (0 == c && 0 == strcmp("pen-gap-extension", options[option_index].name))) {       
          opt->pen_gape = atoi(optarg);
      }
      else if(c == 'H' || (0 == c && 0 == strcmp("vsw-type", options[option_index].name))) {       
          opt->vsw_type = atoi(optarg); 
      }
      else if(c == 'I' || (0 == c && 0 == strcmp("use-seq-equal", options[option_index].name))) {       
          opt->seq_eq = 1;
      }
      else if(c == 'M' || (0 == c && 0 == strcmp("pen-mismatch", options[option_index].name))) {       
          opt->pen_mm = atoi(optarg);
      }                                           
    
      else if(c == 'O' || (0 == c && 0 == strcmp("pen-gap-open", options[option_index].name))) {       
          opt->pen_gapo = atoi(optarg);
      }
      else if(c == 'R' || (0 == c && 0 == strcmp("sam-read-group", options[option_index].name))) {       
          tmap_map_opt_add_tabbed(&opt->sam_rg, &opt->sam_rg_num, optarg);
      }
      else if(c == 'T' || (0 == c && 0 == strcmp("score-thres", options[option_index].name))) {       
          opt->score_thr = atoi(optarg);
      }
      else if(c == 'W' || (0 == c && 0 == strcmp("duplicate-window", options[option_index].name))) {       
          opt->dup_window = atoi(optarg);
      }
      else if(c == 'a' || (0 == c && 0 == strcmp("aln-output-mode", options[option_index].name))) {       
          opt->aln_output_mode = atoi(optarg);
      }
      else if(c == 'f' || (0 == c && 0 == strcmp("fn-fasta", options[option_index].name))) {       
          free(opt->fn_fasta);
          opt->fn_fasta = tmap_strdup(optarg);
      }
      else if(c == 'g' || (0 == c && 0 == strcmp("softclip-type", options[option_index].name))) {       
          opt->softclip_type = atoi(optarg);
      }
      else if(c == 'h' || (0 == c && 0 == strcmp("help", options[option_index].name))) {       
          break;
      }
      else if(c == 'i' || (0 == c && 0 == strcmp("reads-format", options[option_index].name))) {       
          opt->reads_format = tmap_get_reads_file_format_int(optarg);
      }
      else if(c == 'j' || (0 == c && 0 == strcmp("input-bz2", options[option_index].name))) {       
          opt->input_compr = TMAP_FILE_BZ2_COMPRESSION;
          for(i=0;i<opt->fn_reads_num;i++) {
              tmap_get_reads_file_format_from_fn_int(opt->fn_reads[i], &opt->reads_format, &opt->input_compr);
          }
      }
      else if(c == 'k' || (0 == c && 0 == strcmp("shared-memory-key", options[option_index].name))) {       
          opt->shm_key = atoi(optarg);
      }
      else if(c == 'n' || (0 == c && 0 == strcmp("num-threads", options[option_index].name))) {       
          opt->num_threads = atoi(optarg);
          opt->num_threads_autodetected = 0;
      }
      else if(c == 'q' || (0 == c && 0 == strcmp("reads-queue-size", options[option_index].name))) {       
          opt->reads_queue_size = atoi(optarg);
      }
      else if(c == 'r' || (0 == c && 0 == strcmp("fn-reads", options[option_index].name))) {       
          opt->fn_reads_num++;
          opt->fn_reads = tmap_realloc(opt->fn_reads, sizeof(char*) * opt->fn_reads_num, "opt->fn_reads");
          opt->fn_reads[opt->fn_reads_num-1] = tmap_strdup(optarg); 
          tmap_get_reads_file_format_from_fn_int(opt->fn_reads[opt->fn_reads_num-1], &opt->reads_format, &opt->input_compr);
      }
      else if(c == 's' || (0 == c && 0 == strcmp("fn-sam", options[option_index].name))) {       
          free(opt->fn_sam);
          opt->fn_sam = tmap_strdup(optarg);
      }
      else if(c == 'v' || (0 == c && 0 == strcmp("verbose", options[option_index].name))) {       
          tmap_progress_set_verbosity(1);
      }
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
      else if(c == 'x' || (0 == c && 0 == strcmp("sample-reads", options[option_index].name))) {
          opt->sample_reads = atof(optarg); 
      }
#endif
      else if(c == 'w' || (0 == c && 0 == strcmp("band-width", options[option_index].name))) {       
          opt->bw = atoi(optarg);
      }
      else if(c == 'z' || (0 == c && 0 == strcmp("input-gz", options[option_index].name))) {       
          opt->input_compr = TMAP_FILE_GZ_COMPRESSION;
          for(i=0;i<opt->fn_reads_num;i++) {
              tmap_get_reads_file_format_from_fn_int(opt->fn_reads[i], &opt->reads_format, &opt->input_compr);
          }
      }
      // End of global options
      // Flowspace options
      else if(c == 'F' || (0 == c && 0 == strcmp("final-flowspace", options[option_index].name))) {       
          opt->aln_flowspace = 1;
      }
      else if(c == 'N' || (0 == c && 0 == strcmp("use-flowgram", options[option_index].name))) {       
          opt->ignore_flowgram = 1;
      }
      else if(c == 'X' || (0 == c && 0 == strcmp("pen-flow-error", options[option_index].name))) {       
          opt->fscore = atoi(optarg);
      }
      else if(c == 'Y' || (0 == c && 0 == strcmp("sam-flowspace-tags", options[option_index].name))) {       
          opt->sam_flowspace_tags = 1;
      }
      else if(c == 'y' || (0 == c && 0 == strcmp("softclip-key", options[option_index].name))) {       
          opt->softclip_key = 1;
      }
      // End of flowspace options
      // Pairing options
      else if(c == 'Q' || (0 == c && 0 == strcmp("pairing", options[option_index].name))) {
          opt->pairing = atoi(optarg);
          if(1 == opt->pairing) {
              opt->strandedness = 0;
              opt->positioning = 1;
          }
          else if(2 == opt->pairing) {
              opt->strandedness = 1;
              opt->positioning = 0;
          }
      }
      else if(c == 'S' || (0 == c && 0 == strcmp("strandedness", options[option_index].name))) {
          opt->strandedness = atoi(optarg);
      }
      else if(c == 'P' || (0 == c && 0 == strcmp("positioning", options[option_index].name))) {
          opt->positioning = atoi(optarg);
      }
      else if(c == 'b' || (0 == c && 0 == strcmp("ins-size-mean", options[option_index].name))) {
          opt->ins_size_mean = atof(optarg);
      }
      else if(c == 'c' || (0 == c && 0 == strcmp("ins-size-std", options[option_index].name))) {
          opt->ins_size_std = atof(optarg);
      }
      else if(c == 'd' || (0 == c && 0 == strcmp("ins-size-std-max-num", options[option_index].name))) {
          opt->ins_size_std_max_num = atof(optarg);
      }
      else if(c == 'L' || (0 == c && 0 == strcmp("read-rescue", options[option_index].name))) {
          opt->read_rescue = 1;
      }
      else if(c == 'l' || (0 == c && 0 == strcmp("read-rescue-std-num", options[option_index].name))) {
          opt->read_rescue_std_num = atof(optarg);
      }
      else if(c == 'm' || (0 == c && 0 == strcmp("read-rescue-mapq-thr", options[option_index].name))) {
          opt->read_rescue_mapq_thr = atoi(optarg);
      }
      // End of pairing options 
      // End single flag options
      else if(0 != c) {
          tmap_bug();
      }
      // MAP1/MAP2/MAP3/MAPVSW
      else if(0 == strcmp("min-seq-length", options[option_index].name) && (opt->algo_id == TMAP_MAP_ALGO_MAP1 || opt->algo_id == TMAP_MAP_ALGO_MAP2 
                                                                            || opt->algo_id == TMAP_MAP_ALGO_MAP3 || opt->algo_id == TMAP_MAP_ALGO_MAPVSW)) {
          opt->min_seq_len = atoi(optarg);
      }
      else if(0 == strcmp("max-seq-length", options[option_index].name) && (opt->algo_id == TMAP_MAP_ALGO_MAP1 || opt->algo_id == TMAP_MAP_ALGO_MAP2 
                                                                            || opt->algo_id == TMAP_MAP_ALGO_MAP3 || opt->algo_id == TMAP_MAP_ALGO_MAPVSW)) {
          opt->max_seq_len = atoi(optarg);
      }
      // MAP1/MAP3
      else if(0 == strcmp("seed-length", options[option_index].name) && (opt->algo_id == TMAP_MAP_ALGO_MAP1 || opt->algo_id == TMAP_MAP_ALGO_MAP3)) {
          opt->seed_length = atoi(optarg);
          opt->seed_length_set = 1;
      }
      // MAP2/MAP3
      else if(0 == strcmp("max-seed-hits", options[option_index].name) && (opt->algo_id == TMAP_MAP_ALGO_MAP2 || opt->algo_id == TMAP_MAP_ALGO_MAP3)) {
          opt->max_seed_hits = atoi(optarg);
      }
      // MAP3/MAP4
      else if(0 == strcmp("hit-frac", options[option_index].name) && (opt->algo_id == TMAP_MAP_ALGO_MAP3 || opt->algo_id == TMAP_MAP_ALGO_MAP4)) {
          opt->hit_frac = atof(optarg);
      }
      else if(0 == strcmp("seed-step", options[option_index].name) && (opt->algo_id == TMAP_MAP_ALGO_MAP3 || opt->algo_id == TMAP_MAP_ALGO_MAP4)) {
          opt->seed_step = atoi(optarg);
      }
      // MAP1
      else if(0 == strcmp("seed-max-diff", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          opt->seed_max_diff = atoi(optarg);
      }
      else if(0 == strcmp("seed2-length", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          opt->seed2_length = atoi(optarg);
      }
      else if(0 == strcmp("max-diff", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          if(NULL != strstr(optarg, ".")) opt->max_diff = -1, opt->max_diff_fnr = atof(optarg);
          else opt->max_diff = atoi(optarg), opt->max_diff_fnr = -1.0;
      }
      else if(0 == strcmp("max-error-rate", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          opt->max_err_rate = atof(optarg);
      }
      else if(0 == strcmp("max-mismatches", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          if(NULL != strstr(optarg, ".")) opt->max_mm = -1, opt->max_mm_frac = atof(optarg);
          else opt->max_mm = atoi(optarg), opt->max_mm_frac = -1.0;
      }
      else if(0 == strcmp("max-gap-opens", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          if(NULL != strstr(optarg, ".")) opt->max_gapo = -1, opt->max_gapo_frac = atof(optarg);
          else opt->max_gapo = atoi(optarg), opt->max_gapo_frac = -1.0;
      }
      else if(0 == strcmp("max-gap-extensions", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          if(NULL != strstr(optarg, ".")) opt->max_gape = -1, opt->max_gape_frac = atof(optarg);
          else opt->max_gape = atoi(optarg), opt->max_gape_frac = -1.0;
      }
      else if(0 == strcmp("max-cals-deletion", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          opt->max_cals_del = atoi(optarg);
      }
      else if(0 == strcmp("indel-ends-bound", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          opt->indel_ends_bound = atoi(optarg);
      }
      else if(0 == strcmp("max-best-cals", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          opt->max_best_cals = atoi(optarg);
      }
      else if(0 == strcmp("max-nodes", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP1) {
          opt->max_entries = atoi(optarg);
      }
      // MAP2
      else if(0 == strcmp("length-coef", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP2) {
          opt->length_coef = atof(optarg);
      }
      else if(0 == strcmp("max-seed-intv", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP2) {
          opt->max_seed_intv = atoi(optarg);
      }
      else if(0 == strcmp("z-best", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP2) {
          opt->z_best= atoi(optarg);
      }
      else if(0 == strcmp("seeds-rev", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP2) {
          opt->seeds_rev = atoi(optarg);
      }
      else if(0 == strcmp("narrow-rmdup", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP2) {
          opt->narrow_rmdup = 1;
      }
      // MAP 3
      else if(0 == strcmp("hp-diff", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP3) {
          opt->hp_diff = atoi(optarg);
      }
      else if(0 == strcmp("fwd-search", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP3) {
          opt->fwd_search = 1;
      }
      else if(0 == strcmp("skip-seed-frac", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP3) {
          opt->skip_seed_frac = atof(optarg);
      }
      // MAP 4
      else if(0 == strcmp("min-seed-length", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP4) {
          opt->min_seed_length= atoi(optarg);
      }
      else if(0 == strcmp("max-seed-length", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP4) {
          opt->max_seed_length= atoi(optarg);
      }
      else if(0 == strcmp("max-iwidth", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP4) {
          opt->max_iwidth = atoi(optarg);
      }
      else if(0 == strcmp("max-repr", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP4) {
          opt->max_repr = atoi(optarg);
      }
      else if(0 == strcmp("rand-repr", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP4) {
          opt->rand_repr = 1;
      }
      else if(0 == strcmp("use-min", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_MAP4) {
          opt->use_min = 1;
      }
      // STAGE
      else if(0 == strcmp("stage-score-thres", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_STAGE) {
          opt->stage_score_thr = atoi(optarg);
      }
      else if(0 == strcmp("stage-mapq-thres", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_STAGE) {
          opt->stage_mapq_thr = atoi(optarg);
      }
      else if(0 == strcmp("stage-keep-all", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_STAGE) {
          opt->stage_keep_all = 1;
      }
      else if(0 == strcmp("stage-seed-freq-cutoff", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_STAGE) {
          opt->stage_seed_freqc = atof(optarg);
      }
      else if(0 == strcmp("stage-seed-freq-cutoff-group-frac", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_STAGE) {
          opt->stage_seed_freqc_group_frac = atof(optarg);
      }
      else if(0 == strcmp("stage-seed-freq-cutoff-rand-repr", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_STAGE) {
          opt->stage_seed_freqc_rand_repr = atoi(optarg);
      }
      else if(0 == strcmp("stage-seed-freq-cutoff-min-groups", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_STAGE) {
          opt->stage_seed_freqc_min_groups = atoi(optarg);
      }
      else if(0 == strcmp("stage-seed-max-length", options[option_index].name) && opt->algo_id == TMAP_MAP_ALGO_STAGE) {
          opt->stage_seed_max_length = atoi(optarg);
      }
      // MAPALL
      
      else {
          break;
      }
      if(argc == optind) {
          val = 1;
      }
  }
  if(optind < argc) {
      val = 0;
      tmap_file_fprintf(tmap_file_stderr, "non-option command-line-elements:\n");
      while(optind < argc) {
          tmap_file_fprintf(tmap_file_stderr, "%s", argv[optind]);
          i = -1;
          if(opt->algo_id == TMAP_MAP_ALGO_MAPALL) {
              i = tmap_algo_name_to_id(argv[optind]);
              if(0 <= i) {
                  tmap_file_fprintf(tmap_file_stderr, ": recognized an algorithm name, did you forget to include the stage parameter?\n");
              }
          }
          if(i < 0) {
              tmap_file_fprintf(tmap_file_stderr, ": unknown command line option\n");
          }
          optind++;
      }
  }
  free(options);
  free(getopt_format);
  return val;
}

static int32_t
tmap_map_opt_file_check_with_null(char *fn1, char *fn2)
{
  if(NULL == fn1 && NULL == fn2) {
      return 0;
  }
  else if((NULL == fn1 && NULL != fn2)
          || (NULL != fn1 && NULL == fn2)) {
      return 1;
  }
  else if(0 != strcmp(fn1, fn2)) {
      return 1;
  }
  return 0;
}

// check that the global and flowspace options match the algorithm specific global options
void
tmap_map_opt_check_global(tmap_map_opt_t *opt_a, tmap_map_opt_t *opt_b) 
{
    int32_t i;
    // global
    if(0 != tmap_map_opt_file_check_with_null(opt_a->fn_fasta, opt_b->fn_fasta)) {
        tmap_error("option -f was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->fn_reads_num != opt_b->fn_reads_num) {
        tmap_error("option -r was specified outside of the common options", Exit, CommandLineArgument);
    }
    for(i=0;i<opt_a->fn_reads_num;i++) {
        if(0 != tmap_map_opt_file_check_with_null(opt_a->fn_reads[i], opt_b->fn_reads[i])) {
            tmap_error("option -r was specified outside of the common options", Exit, CommandLineArgument);
        }
    }
    if(0 != tmap_map_opt_file_check_with_null(opt_a->fn_sam, opt_b->fn_sam)) {
        tmap_error("option -s was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->reads_format != opt_b->reads_format) {
        tmap_error("option -i was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->score_match != opt_b->score_match) {
        tmap_error("option -A was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->pen_mm != opt_b->pen_mm) {
        tmap_error("option -M was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->pen_gapo != opt_b->pen_gapo) {
        tmap_error("option -O was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->pen_gape != opt_b->pen_gape) {
        tmap_error("option -E was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->bw != opt_b->bw) {
        tmap_error("option -w was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->softclip_type != opt_b->softclip_type) {
        tmap_error("option -g was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->dup_window != opt_b->dup_window) {
        tmap_error("option -W was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->max_seed_band != opt_b->max_seed_band) {
        tmap_error("option -B was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->unroll_banding != opt_b->unroll_banding) {
        tmap_error("option -U was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->score_thr != opt_b->score_thr) {
        tmap_error("option -T was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->reads_queue_size != opt_b->reads_queue_size) {
        tmap_error("option -q was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->num_threads != opt_b->num_threads) {
        tmap_error("option -n was specified outside of the common options", Exit, CommandLineArgument);
    }
    /* NB: "aln_output_mode" or "-a" may be modified by mapall */
    if(opt_a->sam_rg_num != opt_b->sam_rg_num) {
        tmap_error("option -R was specified outside of the common options", Exit, CommandLineArgument);
    }
    for(i=0;i<opt_a->sam_rg_num;i++) {
        if(0 != tmap_map_opt_file_check_with_null(opt_a->sam_rg[i], opt_b->sam_rg[i])) {
            tmap_error("option -R was specified outside of the common options", Exit, CommandLineArgument);
        }
    }
    if(opt_a->bidirectional != opt_b->bidirectional) {
        tmap_error("option -D was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->seq_eq != opt_b->seq_eq) {
        tmap_error("option -I was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->input_compr != opt_b->input_compr) {
        tmap_error("option -j or -z was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->shm_key != opt_b->shm_key) {
        tmap_error("option -k was specified outside of the common options", Exit, CommandLineArgument);
    }
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
    if(opt_a->sample_reads != opt_b->sample_reads) {
        tmap_error("option -x was specified outside of the common options", Exit, CommandLineArgument);
    }
#endif
    if(opt_a->vsw_type != opt_b->vsw_type) {
        tmap_error("option -H was specified outside of the common options", Exit, CommandLineArgument);
    }
    // flowspace
    if(opt_a->fscore != opt_b->fscore) {
        tmap_error("option -X was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->softclip_key != opt_b->softclip_key) {
        tmap_error("option -y was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->sam_flowspace_tags != opt_b->sam_flowspace_tags) {
        tmap_error("option -Y was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->ignore_flowgram != opt_b->ignore_flowgram) {
        tmap_error("option -S was specified outside of the common options", Exit, CommandLineArgument);
    }
    if(opt_a->aln_flowspace != opt_b->aln_flowspace) {
        tmap_error("option -F was specified outside of the common options", Exit, CommandLineArgument);
    }
    // pairing
    if(opt_a->pairing != opt_b->pairing) {
        tmap_error("option -Q was specified outside the common options", Exit, CommandLineArgument);
    }
    if(opt_a->strandedness != opt_b->strandedness) {
        tmap_error("option -S was specified outside the common options", Exit, CommandLineArgument);
    }
    if(opt_a->positioning != opt_b->positioning) {
        tmap_error("option -P was specified outside the common options", Exit, CommandLineArgument);
    }
    if(opt_a->ins_size_mean != opt_b->ins_size_mean) {
        tmap_error("option -b specified outside the common options", Exit, CommandLineArgument);
    }
    if(opt_a->ins_size_std != opt_b->ins_size_std) {
        tmap_error("option -c was specified outside the common options", Exit, CommandLineArgument);
    }
    if(opt_a->ins_size_std_max_num != opt_b->ins_size_std_max_num) {
        tmap_error("option -d was specified outside the common options", Exit, CommandLineArgument);
    }
    if(opt_a->read_rescue != opt_b->read_rescue) {
        tmap_error("option -L was specified outside the common options", Exit, CommandLineArgument);
    }
    if(opt_a->read_rescue_std_num != opt_b->read_rescue_std_num) {
        tmap_error("option -l was specified outside the common options", Exit, CommandLineArgument);
    }
    if(opt_a->read_rescue_mapq_thr != opt_b->read_rescue_mapq_thr) {
        tmap_error("option -m was specified outside the common options", Exit, CommandLineArgument);
    }
}

void
tmap_map_opt_check_stage(tmap_map_opt_t *opt_a, tmap_map_opt_t *opt_b) 
{
  if(opt_a->stage_score_thr != opt_b->stage_score_thr) {
      tmap_error("option --stage-score-thres specified outside of stage options", Exit, CommandLineArgument);
  }
  if(opt_a->stage_mapq_thr != opt_b->stage_mapq_thr) {
      tmap_error("option --stage-mapq-thres specified outside of stage options", Exit, CommandLineArgument);
  }
  if(opt_a->stage_keep_all != opt_b->stage_keep_all) {
      tmap_error("option --stage-keep-all specified outside of stage options", Exit, CommandLineArgument);
  }
  if(opt_a->stage_seed_freqc != opt_b->stage_seed_freqc) {
      tmap_error("option --stage-seed-freqc specified outside of stage options", Exit, CommandLineArgument);
  }
  if(opt_a->stage_seed_freqc_group_frac != opt_b->stage_seed_freqc_group_frac) {
      tmap_error("option --stage-seed-freqc-group-frac specified outside of stage options", Exit, CommandLineArgument);
  }
  if(opt_a->stage_seed_freqc_rand_repr != opt_b->stage_seed_freqc_rand_repr) {
      tmap_error("option --stage-seed-freqc-rand-repr specified outside of stage options", Exit, CommandLineArgument);
  }
  if(opt_a->stage_seed_freqc_min_groups != opt_b->stage_seed_freqc_min_groups) {
      tmap_error("option --stage-seed-freqc-min-groups specified outside of stage options", Exit, CommandLineArgument);
  }
  if(opt_a->stage_seed_max_length != opt_b->stage_seed_max_length) {
      tmap_error("option --stage-score-thres specified outside of stage options", Exit, CommandLineArgument);
  }
}

void
tmap_map_opt_check(tmap_map_opt_t *opt)
{
  //int32_t i;
  // global and flowspace options
  if(NULL == opt->fn_fasta && 0 == opt->shm_key) {
      tmap_error("option -f or option -k must be specified", Exit, CommandLineArgument);
  }
  else if(NULL != opt->fn_fasta && 0 < opt->shm_key) {
      tmap_error("option -f and option -k may not be specified together", Exit, CommandLineArgument);
  }
  if(0 == opt->fn_reads_num && TMAP_READS_FORMAT_UNKNOWN == opt->reads_format) {
      tmap_error("option -r or option -i must be specified", Exit, CommandLineArgument);
  }
  else if(1 < opt->fn_reads_num) {
      if(1 == opt->sam_flowspace_tags) {
          tmap_error("options -1 and -2 cannot be used with -Y", Exit, CommandLineArgument);
      }
      else if(0 != opt->pairing) {
          if(opt->strandedness < 0 || 1 < opt->strandedness) {
              tmap_error("option -S was not specified", Exit, CommandLineArgument);
          }
          else if(opt->positioning < 0 || 1 < opt->positioning) {
              tmap_error("option -P was not specified", Exit, CommandLineArgument);
          }
          else if(opt->ins_size_mean < 0) {
              tmap_error("option -b not specified", Exit, CommandLineArgument);
          }
          else if(opt->ins_size_std < 0) {
              tmap_error("option -c not specified", Exit, CommandLineArgument);
          }
          else if(opt->ins_size_std_max_num < 0) {
              tmap_error("option -d was not specified", Exit, CommandLineArgument);
          }
          else if(1 == opt->read_rescue) {
              if(opt->read_rescue_std_num < 0) {
                  tmap_error("option -l was not specified", Exit, CommandLineArgument);
              }
              tmap_error_cmd_check_int(opt->read_rescue_mapq_thr, 0, 1, "-m");
          }
      }
      // OK
  }
  if(TMAP_READS_FORMAT_UNKNOWN == opt->reads_format) {
      tmap_error("the reads format (-r/-i) was unrecognized", Exit, CommandLineArgument);
  }
  tmap_error_cmd_check_int(opt->score_match, 1, INT32_MAX, "-A");
  tmap_error_cmd_check_int(opt->pen_mm, 1, INT32_MAX, "-M");
  tmap_error_cmd_check_int(opt->pen_gapo, 1, INT32_MAX, "-O");
  tmap_error_cmd_check_int(opt->pen_gape, 1, INT32_MAX, "-E");
  tmap_error_cmd_check_int(opt->fscore, 0, INT32_MAX, "-X");
  tmap_error_cmd_check_int(opt->bw, 0, INT32_MAX, "-w");
  tmap_error_cmd_check_int(opt->softclip_type, 0, 3, "-g");
  tmap_error_cmd_check_int(opt->softclip_key, 0, 0, "-y");

  tmap_error_cmd_check_int(opt->dup_window, -1, INT32_MAX, "-W");
  tmap_error_cmd_check_int(opt->max_seed_band, 1, INT32_MAX, "-B");
  tmap_error_cmd_check_int(opt->unroll_banding, 0, 1, "-U");
  tmap_error_cmd_check_int(opt->score_thr, INT32_MIN, INT32_MAX, "-T");
  if(-1 != opt->reads_queue_size) tmap_error_cmd_check_int(opt->reads_queue_size, 1, INT32_MAX, "-q");
  tmap_error_cmd_check_int(opt->num_threads, 1, INT32_MAX, "-n");
  tmap_error_cmd_check_int(opt->aln_output_mode, 0, 3, "-a");
  tmap_error_cmd_check_int(opt->bidirectional, 0, 1, "-D");
  tmap_error_cmd_check_int(opt->seq_eq, 0, 1, "-I");
  tmap_error_cmd_check_int(opt->ignore_rg_sam_tags, 0, 1, "-C");
  /*
  if(0 == opt->ignore_rg_sam_tags && NULL != opt->sam_rg) {
      tmap_error("Must use -C with -R", Exit, CommandLineArgument);
  }
  */
  if(-1 != opt->min_seq_len) tmap_error_cmd_check_int(opt->min_seq_len, 1, INT32_MAX, "--min-seq-length");
  if(-1 != opt->max_seq_len) tmap_error_cmd_check_int(opt->max_seq_len, 1, INT32_MAX, "--max-seq-length");
  if(-1 != opt->min_seq_len && -1 != opt->max_seq_len && opt->max_seq_len < opt->min_seq_len) {
      tmap_error("The minimum sequence length must be less than the maximum sequence length (--min-seq-length and --max-seq-length)", Exit, CommandLineArgument);
  }
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
  tmap_error_cmd_check_int(opt->sample_reads, 0, 1, "-x");
#endif
  tmap_error_cmd_check_int(opt->vsw_type, 1, 10, "-H");
  // Warn users
  switch(opt->vsw_type) {
    case 1:
    case 4:
    case 6:
      break;
    default:
      tmap_error("the option -H value has not been extensively tested; proceed with caution", Warn, CommandLineArgument);
      break;
  }

  switch(opt->algo_id) {
    case TMAP_MAP_ALGO_MAP1: // map1 options
      if(-1 != opt->seed_length) tmap_error_cmd_check_int(opt->seed_length, 1, INT32_MAX, "--seed-length");
      tmap_error_cmd_check_int(opt->seed_max_diff, 0, INT32_MAX, "--seed-max-diff");
      if(-1 != opt->seed2_length) tmap_error_cmd_check_int(opt->seed2_length, 1, INT32_MAX, "--seed2-length");
      if(-1 != opt->seed_length && -1 != opt->seed2_length) {
          tmap_error_cmd_check_int(opt->seed_length, 1, opt->seed2_length, "The secondary seed length (--seed2-length) must be greater than the primary seed length (--seed-length)");
      }
      tmap_error_cmd_check_int((opt->max_diff_fnr < 0) ? opt->max_diff: (int32_t)opt->max_diff_fnr, 0, INT32_MAX, "--max-diff");
      tmap_error_cmd_check_int((int32_t)opt->max_err_rate, 0, INT32_MAX, "--max-error-rate");
      // this will take care of the case where they are both < 0
      tmap_error_cmd_check_int((opt->max_mm_frac < 0) ? opt->max_mm : (int32_t)opt->max_mm_frac, 0, INT32_MAX, "--max-mismatches");
      // this will take care of the case where they are both < 0
      tmap_error_cmd_check_int((opt->max_gapo_frac < 0) ? opt->max_gapo : (int32_t)opt->max_gapo_frac, 0, INT32_MAX, "--max-gap-opens");
      // this will take care of the case where they are both < 0
      tmap_error_cmd_check_int((opt->max_gape_frac < 0) ? opt->max_gape : (int32_t)opt->max_gape_frac, 0, INT32_MAX, "--max-gap-extensions");
      tmap_error_cmd_check_int(opt->max_cals_del, 1, INT32_MAX, "--max-cals-deletion");
      tmap_error_cmd_check_int(opt->indel_ends_bound, 0, INT32_MAX, "--indels-ends-bound");
      tmap_error_cmd_check_int(opt->max_best_cals, 0, INT32_MAX, "--max-best-cals");
      tmap_error_cmd_check_int(opt->max_entries, 1, INT32_MAX, "--max-nodes");
      break;
    case TMAP_MAP_ALGO_MAP2:
      //tmap_error_cmd_check_int(opt->mask_level, 0, 1, "-m");
      tmap_error_cmd_check_int(opt->length_coef, 0, INT32_MAX, "--length-coef");
      tmap_error_cmd_check_int(opt->max_seed_intv, 0, INT32_MAX, "--max-seed-intv");
      tmap_error_cmd_check_int(opt->z_best, 1, INT32_MAX, "--z-best");
      tmap_error_cmd_check_int(opt->seeds_rev, 0, INT32_MAX, "--seeds-rev");
      tmap_error_cmd_check_int(opt->narrow_rmdup, 0, 1, "--narrow-rmdup");
      tmap_error_cmd_check_int(opt->max_seed_hits, 1, INT32_MAX, "--max-seed-hits");
      break;
    case TMAP_MAP_ALGO_MAP3:
      if(-1 != opt->seed_length) tmap_error_cmd_check_int(opt->seed_length, 1, INT32_MAX, "--seed-length");
      tmap_error_cmd_check_int(opt->max_seed_hits, 1, INT32_MAX, "--max-seed-hits");
      tmap_error_cmd_check_int(opt->hp_diff, 0, INT32_MAX, "--hp-diff");
      tmap_error_cmd_check_int(opt->hit_frac, 0, 1, "--hit-frac");
      if(-1 != opt->seed_step) tmap_error_cmd_check_int(opt->seed_step, 1, INT32_MAX, "--seed-step");
      tmap_error_cmd_check_int(opt->skip_seed_frac, 0, 1, "--skip-seed-frac");
      break;
    case TMAP_MAP_ALGO_MAP4:
      if(-1 != opt->min_seed_length) {
          tmap_error_cmd_check_int(opt->min_seed_length, 1, INT32_MAX, "--min-seed-length");
          if(opt->max_seed_length < opt->min_seed_length) {
              tmap_error("--max-seed-length is less than --min-seed-length", Exit, CommandLineArgument);
          }
      }
      if(-1 != opt->max_seed_length) tmap_error_cmd_check_int(opt->max_seed_length, 1, INT32_MAX, "--max-seed-length");
      tmap_error_cmd_check_int(opt->hit_frac, 0, 1, "--hit-frac");
      if(-1 != opt->seed_step) tmap_error_cmd_check_int(opt->seed_step, 1, INT32_MAX, "--seed-step");
      tmap_error_cmd_check_int(opt->max_iwidth, 0, INT32_MAX, "--max-iwidth");
      tmap_error_cmd_check_int(opt->max_repr, 0, INT32_MAX, "--max-repr");
      tmap_error_cmd_check_int(opt->rand_repr, 0, 1, "--rand-repr");
      tmap_error_cmd_check_int(opt->use_min, 0, 1, "--use-min");
      break;
    case TMAP_MAP_ALGO_MAPALL:
      if(0 == opt->num_sub_opts) {
          tmap_error("no stages/algorithms given", Exit, CommandLineArgument);
      }
      break;
    case TMAP_MAP_ALGO_STAGE:
      // stage options
      tmap_error_cmd_check_int(opt->stage_score_thr, INT32_MIN, INT32_MAX, "--stage-score-thres");
      tmap_error_cmd_check_int(opt->stage_mapq_thr, 0, 255, "--stage-mapq-thres");
      tmap_error_cmd_check_int(opt->stage_keep_all, 0, 1, "--stage-keep-all");
      tmap_error_cmd_check_int(opt->stage_seed_freqc, 0.0, 1.0, "--stage-seed-freq-cutoff");
      tmap_error_cmd_check_int(opt->stage_seed_freqc_group_frac, 0.0, 1.0, "--stage-seed-freq-cutoff-group-frac");
      tmap_error_cmd_check_int(opt->stage_seed_freqc_rand_repr, 0, INT32_MAX, "--stage-seed-freq-cutoff-rand-repr");
      tmap_error_cmd_check_int(opt->stage_seed_freqc_min_groups, 0, INT32_MAX, "--stage-seed-freq-cutoff-min-groups");
      if(-1 != opt->stage_seed_max_length) tmap_error_cmd_check_int(opt->stage_seed_max_length, 1, INT32_MAX, "--stage-max-seed-length");
      break;
    default:
      break;
  }
  /*
  // check sub-options
  for(i=0;i<opt->num_sub_opts;i++) {
      // check mapping algorithm specific options
      tmap_map_opt_check(opt->sub_opts[i]);

      // check that common values match other opt values
      tmap_map_opt_check_common(opt, opt->sub_opts[i]);
  }
  */
}

void
tmap_map_opt_copy_global(tmap_map_opt_t *opt_dest, tmap_map_opt_t *opt_src)
{
    int i;

    // global options
    opt_dest->fn_fasta = tmap_strdup(opt_src->fn_fasta);
    opt_dest->fn_reads_num = opt_src->fn_reads_num;
    opt_dest->fn_reads = tmap_malloc(sizeof(char*)*opt_dest->fn_reads_num, "opt_dest->fn_reads");
    for(i=0;i<opt_dest->fn_reads_num;i++) {
        opt_dest->fn_reads[i] = tmap_strdup(opt_src->fn_reads[i]);
    }
    opt_dest->reads_format = opt_src->reads_format;
    opt_dest->fn_sam = tmap_strdup(opt_src->fn_sam);
    opt_dest->score_match = opt_src->score_match;
    opt_dest->pen_mm = opt_src->pen_mm;
    opt_dest->pen_gapo = opt_src->pen_gapo;
    opt_dest->pen_gape = opt_src->pen_gape;
    opt_dest->bw = opt_src->bw;
    opt_dest->softclip_type = opt_src->softclip_type;
    opt_dest->dup_window = opt_src->dup_window;
    opt_dest->max_seed_band = opt_src->max_seed_band;
    opt_dest->unroll_banding = opt_src->unroll_banding;
    opt_dest->score_thr = opt_src->score_thr;
    opt_dest->reads_queue_size = opt_src->reads_queue_size;
    opt_dest->num_threads = opt_src->num_threads;
    if(TMAP_MAP_ALGO_STAGE == opt_dest->algo_id) {
        opt_dest->aln_output_mode = opt_src->aln_output_mode;
    }
    else {
        opt_dest->aln_output_mode = TMAP_MAP_OPT_ALN_MODE_ALL;
    }
    opt_dest->sam_rg_num = opt_src->sam_rg_num; 
    opt_dest->sam_rg = tmap_malloc(sizeof(char*)*opt_dest->sam_rg_num, "opt_dest->sam_rg");
    for(i=0;i<opt_dest->sam_rg_num;i++) {
        opt_dest->sam_rg[i] = tmap_strdup(opt_src->sam_rg[i]);
    }
    opt_dest->bidirectional = opt_src->bidirectional;
    opt_dest->seq_eq = opt_src->seq_eq;
    opt_dest->ignore_rg_sam_tags = opt_src->ignore_rg_sam_tags;
    opt_dest->input_compr = opt_src->input_compr;
    opt_dest->shm_key = opt_src->shm_key;
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
    opt_dest->sample_reads = opt_src->sample_reads;
#endif
    opt_dest->vsw_type = opt_src->vsw_type;
    
    // flowspace options
    opt_dest->fscore = opt_src->fscore;
    opt_dest->softclip_key = opt_src->softclip_key;
    opt_dest->sam_flowspace_tags = opt_src->sam_flowspace_tags;
    opt_dest->ignore_flowgram = opt_src->ignore_flowgram;
    opt_dest->aln_flowspace = opt_src->aln_flowspace;

    // pairing
    opt_dest->pairing = opt_src->pairing;
    opt_dest->strandedness = opt_src->strandedness;
    opt_dest->positioning = opt_src->positioning;
    opt_dest->ins_size_mean = opt_src->ins_size_mean;
    opt_dest->ins_size_std = opt_src->ins_size_std;
    opt_dest->ins_size_std_max_num = opt_src->ins_size_std_max_num;
    opt_dest->read_rescue = opt_src->read_rescue;
    opt_dest->read_rescue_std_num = opt_src->read_rescue_std_num;
    opt_dest->read_rescue_mapq_thr = opt_src->read_rescue_mapq_thr;
}

void
tmap_map_opt_copy_stage(tmap_map_opt_t *opt_dest, tmap_map_opt_t *opt_src)
{
  opt_dest->stage_score_thr = opt_src->stage_score_thr;
  opt_dest->stage_mapq_thr = opt_src->stage_mapq_thr;
  opt_dest->stage_keep_all = opt_src->stage_keep_all;
  opt_dest->stage_seed_freqc = opt_src->stage_seed_freqc;
  opt_dest->stage_seed_freqc_group_frac = opt_src->stage_seed_freqc_group_frac;
  opt_dest->stage_seed_freqc_rand_repr = opt_src->stage_seed_freqc_rand_repr;
  opt_dest->stage_seed_freqc_min_groups = opt_src->stage_seed_freqc_min_groups;
  opt_dest->stage_seed_max_length = opt_src->stage_seed_max_length;
}

void
tmap_map_opt_print(tmap_map_opt_t *opt)
{
  int32_t i;
  fprintf(stderr, "algo_id=%d\n", opt->algo_id);
  fprintf(stderr, "fn_fasta=%s\n", opt->fn_fasta);
  for(i=0;i<opt->fn_reads_num;i++) {
      if(0 < i) fprintf(stderr, ",");
      fprintf(stderr, "fn_reads=%s", opt->fn_reads[i]);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "fn_sam=%s\n", opt->fn_sam);
  fprintf(stderr, "reads_format=%d\n", opt->reads_format);
  fprintf(stderr, "reads_format=%d\n", opt->reads_format);
  fprintf(stderr, "score_match=%d\n", opt->score_match);
  fprintf(stderr, "pen_mm=%d\n", opt->pen_mm);
  fprintf(stderr, "pen_gapo=%d\n", opt->pen_gapo);
  fprintf(stderr, "pen_gape=%d\n", opt->pen_gape);
  fprintf(stderr, "fscore=%d\n", opt->fscore);
  fprintf(stderr, "bw=%d\n", opt->bw);
  fprintf(stderr, "softclip_type=%d\n", opt->softclip_type);
  fprintf(stderr, "softclip_key=%d\n", opt->softclip_key);
  fprintf(stderr, "dup_window=%d\n", opt->dup_window);
  fprintf(stderr, "max_seed_band=%d\n", opt->max_seed_band);
  fprintf(stderr, "unroll_banding=%d\n", opt->unroll_banding);
  fprintf(stderr, "score_thr=%d\n", opt->score_thr);
  fprintf(stderr, "reads_queue_size=%d\n", opt->reads_queue_size);
  fprintf(stderr, "num_threads=%d\n", opt->num_threads);
  fprintf(stderr, "aln_output_mode=%d\n", opt->aln_output_mode);
  for(i=0;i<opt->sam_rg_num;i++) {
      if(0 < i) fprintf(stderr, ",");
      fprintf(stderr, "sam_rg=%s", opt->sam_rg[i]);
  }
  fprintf(stderr, "bidirectional=%d\n", opt->bidirectional);
  fprintf(stderr, "seq_eq=%d\n", opt->seq_eq);
  fprintf(stderr, "ignore_rg_sam_tags=%d\n", opt->ignore_rg_sam_tags);
  fprintf(stderr, "sam_flowspace_tags=%d\n", opt->sam_flowspace_tags);
  fprintf(stderr, "ignore_flowgram=%d\n", opt->ignore_flowgram);
  fprintf(stderr, "aln_flowspace=%d\n", opt->aln_flowspace);
  fprintf(stderr, "input_compr=%d\n", opt->input_compr);
  fprintf(stderr, "shm_key=%d\n", (int)opt->shm_key);
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
  fprintf(stderr, "sample_reads=%lf\n", opt->sample_reads);
#endif
  fprintf(stderr, "vsw_type=%d\n", opt->vsw_type);
  fprintf(stderr, "min_seq_len=%d\n", opt->min_seq_len);
  fprintf(stderr, "max_seq_len=%d\n", opt->max_seq_len);
  fprintf(stderr, "seed_length=%d\n", opt->seed_length);
  fprintf(stderr, "seed_length_set=%d\n", opt->seed_length_set);
  fprintf(stderr, "seed_max_diff=%d\n", opt->seed_max_diff);
  fprintf(stderr, "seed2_length=%d\n", opt->seed2_length);
  fprintf(stderr, "max_diff=%d\n", opt->max_diff);
  fprintf(stderr, "max_diff_fnr=%lf\n", opt->max_diff_fnr);
  fprintf(stderr, "max_err_rate=%lf\n", opt->max_err_rate);
  fprintf(stderr, "max_mm=%d\n", opt->max_mm);
  fprintf(stderr, "max_mm_frac=%lf\n", opt->max_mm_frac);
  fprintf(stderr, "max_gapo=%d\n", opt->max_gapo);
  fprintf(stderr, "max_gapo_frac=%lf\n", opt->max_gapo_frac);
  fprintf(stderr, "max_gape=%d\n", opt->max_gape);
  fprintf(stderr, "max_gape_frac=%lf\n", opt->max_gape_frac);
  fprintf(stderr, "max_cals_del=%d\n", opt->max_cals_del);
  fprintf(stderr, "indel_ends_bound=%d\n", opt->indel_ends_bound);
  fprintf(stderr, "max_best_cals=%d\n", opt->max_best_cals);
  fprintf(stderr, "max_entries=%d\n", opt->max_entries);
  fprintf(stderr, "length_coef=%lf\n", opt->length_coef);
  fprintf(stderr, "max_seed_intv=%d\n", opt->max_seed_intv);
  fprintf(stderr, "z_best=%d\n", opt->z_best);
  fprintf(stderr, "seeds_rev=%d\n", opt->seeds_rev);
  fprintf(stderr, "narrow_rmdup=%d\n", opt->narrow_rmdup);
  fprintf(stderr, "max_seed_hits=%d\n", opt->max_seed_hits);
  fprintf(stderr, "hp_diff=%d\n", opt->hp_diff);
  fprintf(stderr, "hit_frac=%lf\n", opt->hit_frac);
  fprintf(stderr, "seed_step=%d\n", opt->seed_step);
  fprintf(stderr, "fwd_search=%d\n", opt->fwd_search);
  fprintf(stderr, "skip_seed_frac=%lf\n", opt->skip_seed_frac);
  fprintf(stderr, "stage_score_thr=%d\n", opt->stage_score_thr);
  fprintf(stderr, "stage_mapq_thr=%d\n", opt->stage_mapq_thr);
  fprintf(stderr, "stage_keep_all=%d\n", opt->stage_keep_all);
  fprintf(stderr, "stage_seed_freqc=%.2f\n", opt->stage_seed_freqc);
  fprintf(stderr, "stage_seed_freqc_group_frac=%.2f\n", opt->stage_seed_freqc_group_frac);
  fprintf(stderr, "stage_seed_freqc_rand_repr=%d\n", opt->stage_seed_freqc_rand_repr);
  fprintf(stderr, "stage_seed_freqc_min_groups=%d\n", opt->stage_seed_freqc_min_groups);
  fprintf(stderr, "stage_seed_max_length=%d\n", opt->stage_seed_max_length);
}
