/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <math.h>
#include <config.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif
#include <unistd.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sam_print.h"
#include "../util/tmap_sort.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt_gen.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_bwt_match.h"
#include "../index/tmap_sa.h"
#include "../io/tmap_seq_io.h"
#include "../server/tmap_shm.h"
#include "../sw/tmap_sw.h"
#include "tmap_map_util.h"
#include "tmap_map_driver.h"
#include "tmap_map1_aux.h"
#include "tmap_map1.h"

static int32_t g_log_n[256];

// sort by max-score
#define __tmap_map1_sam_sort_score_lt(a, b) ((a).score > (b).score)
TMAP_SORT_INIT(tmap_map1_sam_sort_score, tmap_map_sam_t, __tmap_map1_sam_sort_score_lt)

int32_t
tmap_map1_cal_maxdiff(int32_t l, double err, double thres)
{
  double elambda = exp(-l * err);
  double sum, y = 1.0;
  int32_t k, x = 1;
  for(k=1, sum=elambda;k<1000;k++) {
      y *= l * err;
      x *= k;
      sum += elambda * y / x;
      if (1.0 - sum < thres) return k;
  }
  // default
  return 2;
}

void
tmap_map1_print_max_diff(tmap_map_opt_t *opt, int32_t stage)
{
  int32_t i, k, l;

  // initialize
  for(i=0;i<=TMAP_MAP_UTIL_MAX_DIFF_READ_LENGTH;i++) {
      opt->max_diff_table[i] = 0;
  }

  if(opt->max_diff < 0) {
      if(0 < stage) tmap_progress_print("calculating maximum differences in map1 for stage %d", stage);
      else tmap_progress_print("calculating maximum differences in map1");

      for(i = 17, k = 0;i <= TMAP_MAP_UTIL_MAX_DIFF_READ_LENGTH;i++) {
          l = tmap_map1_cal_maxdiff(i, opt->max_err_rate, opt->max_diff_fnr);
          if(l != k ) {
              tmap_progress_print("%dbp reads will have at most %d differences", i, l);
          }
          opt->max_diff_table[i] = l;
          k = l;
      }
  }
  else {
      for(i=0;i <= TMAP_MAP_UTIL_MAX_DIFF_READ_LENGTH;i++) {
          opt->max_diff_table[i] = opt->max_diff;
      }
  }
}


static inline void
tmap_map1_set_g_log_n()
{
  int32_t i;
  for(i=1;i<256;i++) {
      g_log_n[i] = (int32_t)(4.343 * log(i) + 0.5);
  }
}

static inline uint8_t
tmap_map1_sam_mapq(int32_t num_best_sa, int32_t num_all_sa, int32_t max_mm, int32_t num_mm)
{
  int32_t n;

  if(1 < num_best_sa) {
      return 0; // multiple best hits
  }
  else if(max_mm == num_mm) {
      return 25; // maximum possible score
  }

  n = (num_all_sa - num_best_sa);
  if(0 == n) {
      return 37; // no second best hits
  }
  else if(255 < n) {
      n = 255;
  }

  // use MAQ-like mapping qualities
  return (23 < g_log_n[n])? 0 : 23 - g_log_n[n];
}

static inline void
tmap_map1_sams_mapq(tmap_map_sams_t *sams, tmap_map_opt_t *opt)
{
  int32_t i;
  int32_t num_best_sa, num_best, num_all_sa;

  if(0 == sams->n) {
      return;
  }

  // sort by decreasing score
  tmap_sort_introsort(tmap_map1_sam_sort_score, sams->n, sams->sams);

  //Note: assumes that the alignments are sorted by decreasing score
  num_best = num_best_sa = num_all_sa = 0;
  for(i=0;i<sams->n;i++) {
      if(0 < i && sams->sams[i-1].score < sams->sams[i].score) { // check assumption
          tmap_error("bug encountered", Exit, OutOfRange);
      }
      if(sams->sams[i].score < sams->sams[0].score) {
          break;
      }
      num_best++;
      num_best_sa++;
  }
  num_all_sa += sams->n;
  for(i=0;i<num_best;i++) {
      sams->sams[i].mapq = tmap_map1_sam_mapq(num_best_sa, num_all_sa, opt->max_mm, sams->sams[i].aux.map1_aux->n_mm);
  }
  for(i=num_best;i<sams->n;i++) {
      sams->sams[i].mapq = 0;
  }
}

// thread data
typedef struct {
    tmap_bwt_match_width_t *width[2];
    tmap_bwt_match_width_t *seed_width[2];
    int32_t width_length;
    int32_t max_mm;
    int32_t max_gapo;
    int32_t max_gape;
    tmap_map1_aux_stack_t *stack;
} tmap_map1_thread_data_t;

int32_t
tmap_map1_init(tmap_refseq_t *refseq, tmap_map_opt_t *opt)
{
  tmap_map1_print_max_diff(opt, opt->algo_stage);
  opt->score_thr *= opt->score_match;
  // for calculating mapping qualities
  tmap_map1_set_g_log_n();

  return 0;
}

int32_t 
tmap_map1_thread_init(void **data, tmap_map_opt_t *opt)
{
  tmap_map1_thread_data_t *d = NULL;
  d = tmap_calloc(1, sizeof(tmap_map1_thread_data_t), "d");

  d->width[0] = d->width[1] = NULL;
  d->width_length = 0;
  d->stack = NULL;

  d->seed_width[0] = tmap_calloc(1+opt->seed_length, sizeof(tmap_bwt_match_width_t), "seed_width[0]");
  d->seed_width[1] = tmap_calloc(1+opt->seed_length, sizeof(tmap_bwt_match_width_t), "seed_width[1]");

  d->stack = tmap_map1_aux_stack_init();

  // remember to round up
  d->max_mm = (opt->max_mm < 0) ? (int)(0.99 + opt->max_mm_frac * opt->seed2_length) : opt->max_mm; 
  d->max_gapo = (opt->max_gapo < 0) ? (int)(0.99 + opt->max_gapo_frac * opt->seed2_length) : opt->max_gapo; 
  d->max_gape = (opt->max_gape < 0) ? (int)(0.99 + opt->max_gape_frac * opt->seed2_length) : opt->max_gape; 

  (*data) = (void*)d;

  return 0;
}

tmap_map_sams_t*
tmap_map1_thread_map(void **data, tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2], tmap_map_opt_t *opt)
{
  tmap_map1_thread_data_t *d = (tmap_map1_thread_data_t*)(*data);
  int32_t seed2_len = 0, seq_len = 0;;
  tmap_map_opt_t opt_local = (*opt); // copy over values
  tmap_seq_t *seqs[2]={NULL, NULL}, *orig_seq=NULL;
  tmap_string_t *bases[2]={NULL, NULL};
  tmap_map_sams_t *sams = NULL;

  orig_seq = seq;

  if((0 < opt->min_seq_len && tmap_seq_get_bases(orig_seq)->l < opt->min_seq_len)
     || (0 < opt->max_seq_len && opt->max_seq_len < tmap_seq_get_bases(orig_seq)->l)) {
      // go to the next loop
      return tmap_map_sams_init();
  }

  seq_len = tmap_seq_get_bases(orig_seq)->l;

  // not enough bases, ignore
  if(0 < opt->seed_length && seq_len < opt->seed_length){
      return tmap_map_sams_init();
  }

  // clone the sequence 
  seqs[0] = tmap_seq_clone(orig_seq);
  seqs[1] = tmap_seq_clone(orig_seq);

  // Adjust for SFF
  tmap_seq_remove_key_sequence(seqs[0]);
  tmap_seq_remove_key_sequence(seqs[1]);

  tmap_seq_reverse(seqs[0]); // reverse
  tmap_seq_reverse_compliment(seqs[1]); // reverse compliment

  // convert to integers
  tmap_seq_to_int(seqs[0]);
  tmap_seq_to_int(seqs[1]);

  // get bases
  bases[0] = tmap_seq_get_bases(seqs[0]);
  bases[1] = tmap_seq_get_bases(seqs[1]);

  if(opt->seed2_length < 0 || bases[0]->l < opt->seed2_length) {
      seed2_len = seq_len;
      // remember to round up
      opt_local.max_mm = (opt->max_mm < 0) ? (int)(0.99 + opt->max_mm_frac * seed2_len) : opt->max_mm; 
      opt_local.max_gapo = (opt->max_gapo < 0) ? (int)(0.99 + opt->max_gapo_frac * seed2_len) : opt->max_gapo; 
      opt_local.max_gape = (opt->max_gape < 0) ? (int)(0.99 + opt->max_gape_frac * seed2_len) : opt->max_gape; 
  }
  else {
      seed2_len = opt->seed2_length;
      opt_local.max_mm = d->max_mm;
      opt_local.max_gapo = d->max_gapo;
      opt_local.max_gape = d->max_gape;
  }

  // primary width
  if(d->width_length < seq_len) {
      d->width_length = seq_len;
      d->width[0] = tmap_realloc(d->width[0], (1+d->width_length) * sizeof(tmap_bwt_match_width_t), "d->width[0]");
      d->width[1] = tmap_realloc(d->width[1], (1+d->width_length) * sizeof(tmap_bwt_match_width_t), "d->width[1]");
      memset(d->width[0], 0, (1+d->width_length) * sizeof(tmap_bwt_match_width_t));
      memset(d->width[1], 0, (1+d->width_length) * sizeof(tmap_bwt_match_width_t));
  }
  tmap_bwt_match_cal_width_reverse(bwt[0], seq_len, bases[0]->s, d->width[0]);
  tmap_bwt_match_cal_width_reverse(bwt[1], seq_len, bases[1]->s, d->width[1]);

  // seed width
  if(0 < opt->seed_length) {
      tmap_bwt_match_cal_width_reverse(bwt[0], opt->seed_length, bases[0]->s + (seq_len - opt->seed_length), d->seed_width[0]);
      tmap_bwt_match_cal_width_reverse(bwt[1], opt->seed_length, bases[1]->s + (seq_len - opt->seed_length), d->seed_width[1]);
  }

  // map
  sams = tmap_map1_aux_core(seqs, refseq, bwt, sa, d->width, (0 < opt_local.seed_length) ? d->seed_width : NULL, &opt_local, d->stack, seed2_len);

  if(-1 == opt->algo_stage) { // not part of mapall
      // remove duplicates
      tmap_map_util_remove_duplicates(sams, opt->dup_window);

      // mapping quality
      tmap_map1_sams_mapq(sams, opt);

      // filter alignments
      tmap_map_sams_filter(sams, opt->aln_output_mode);

      // re-align the alignments in flow-space
      if(NULL != opt->flow_order) {
          tmap_map_util_fsw(seq,
                            (1 == opt->flow_order_use_sff) ? NULL : opt->flow_order_int,
                            (1 == opt->flow_order_use_sff) ? 0 : strlen(opt->flow_order),
                            sams, refseq, 
                            opt->bw, opt->softclip_type, opt->score_thr,
                            opt->score_match, opt->pen_mm, opt->pen_gapo,
                            opt->pen_gape, opt->fscore);
      }
  }

  // destroy
  tmap_seq_destroy(seqs[0]);
  tmap_seq_destroy(seqs[1]);

  return sams;
}

int32_t
tmap_map1_thread_cleanup(void **data, tmap_map_opt_t *opt)
{
  tmap_map1_thread_data_t *d = (tmap_map1_thread_data_t*)(*data);

  tmap_map1_aux_stack_destroy(d->stack);
  free(d->seed_width[0]);
  free(d->seed_width[1]);
  free(d->width[0]);
  free(d->width[1]);
  free(d);
  (*data) = NULL;
  return 0;
}

static void 
tmap_map1_core(tmap_map_opt_t *opt)
{
  // run the driver
  tmap_map_driver_core(tmap_map1_init, 
                       tmap_map1_thread_init, 
                       tmap_map1_thread_map, 
                       tmap_map1_thread_cleanup,
                       opt);
}

int 
tmap_map1_main(int argc, char *argv[])
{
  tmap_map_opt_t *opt = NULL;

  // init opt
  opt = tmap_map_opt_init(TMAP_MAP_ALGO_MAP1);

  // get options
  if(1 != tmap_map_opt_parse(argc, argv, opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return tmap_map_opt_usage(opt);
  }
  else { 
      // check command line arguments
      tmap_map_opt_check(opt);
  }

  // run map1
  tmap_map1_core(opt);

  // destroy opt
  tmap_map_opt_destroy(opt);

  tmap_progress_print2("terminating successfully");

  return 0;
}
