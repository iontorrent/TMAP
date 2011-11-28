/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <math.h>
#include <config.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif
#include <unistd.h>
#include "../../util/tmap_error.h"
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_definitions.h"
#include "../../util/tmap_progress.h"
#include "../../util/tmap_sam_print.h"
#include "../../util/tmap_sort.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../../index/tmap_bwt_gen.h"
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_bwt_match.h"
#include "../../index/tmap_sa.h"
#include "../../index/tmap_index.h"
#include "../../io/tmap_seq_io.h"
#include "../../server/tmap_shm.h"
#include "../../sw/tmap_sw.h"
#include "../util/tmap_map_stats.h"
#include "../util/tmap_map_util.h"
#include "../tmap_map_driver.h"
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
  for(i=0;i<=TMAP_MAP_OPT_MAX_DIFF_READ_LENGTH;i++) {
      opt->max_diff_table[i] = 0;
  }

  if(opt->max_diff < 0) {
      if(0 < stage) tmap_progress_print("calculating maximum differences in map1 for stage %d", stage);
      else tmap_progress_print("calculating maximum differences in map1");

      for(i = 17, k = 0;i <= TMAP_MAP_OPT_MAX_DIFF_READ_LENGTH;i++) {
          l = tmap_map1_cal_maxdiff(i, opt->max_err_rate, opt->max_diff_fnr);
          if(l != k ) {
              tmap_progress_print("%dbp reads will have at most %d differences", i, l);
          }
          opt->max_diff_table[i] = l;
          k = l;
      }
  }
  else {
      for(i=0;i <= TMAP_MAP_OPT_MAX_DIFF_READ_LENGTH;i++) {
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

static int32_t
tmap_map1_mapq(tmap_map_sams_t *sams, int32_t seq_len, tmap_map_opt_t *opt)
{
  int32_t i;
  int32_t num_best_sa, num_best;

  if(0 == sams->n) {
      return 0;
  }

  // sort by decreasing score
  tmap_sort_introsort(tmap_map1_sam_sort_score, sams->n, sams->sams);

  //Note: assumes that the alignments are sorted by decreasing score
  num_best = num_best_sa = 0;
  for(i=0;i<sams->n;i++) {
      if(0 < i && sams->sams[i-1].score < sams->sams[i].score) { // check assumption
          tmap_bug();
      }
      if(sams->sams[i].score < sams->sams[0].score) {
          break;
      }
      num_best++;
      num_best_sa++;
  }
  for(i=0;i<num_best;i++) {
      sams->sams[i].mapq = tmap_map1_sam_mapq(num_best_sa, 
                                              sams->sams[i].aux.map1_aux->num_all_sa, 
                                              opt->max_mm, 
                                              sams->sams[i].aux.map1_aux->n_mm);
  }
  for(i=num_best;i<sams->n;i++) {
      sams->sams[i].mapq = 0;
  }

  return 0;
}

// thread data
typedef struct {
    tmap_bwt_match_width_t *width;
    tmap_bwt_match_width_t *seed_width;
    int32_t width_length;
    int32_t max_mm;
    int32_t max_gapo;
    int32_t max_gape;
    tmap_map1_aux_stack_t *stack;
} tmap_map1_thread_data_t;

int32_t
tmap_map1_init(void **data, tmap_refseq_t *refseq, tmap_map_opt_t *opt)
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

  d->width = NULL;
  d->width_length = 0;
  d->stack = NULL;

  d->seed_width = tmap_calloc(1+opt->seed_length, sizeof(tmap_bwt_match_width_t), "seed_width");

  d->stack = tmap_map1_aux_stack_init();

  // remember to round up
  d->max_mm = (opt->max_mm < 0) ? (int)(0.99 + opt->max_mm_frac * opt->seed2_length) : opt->max_mm; 
  d->max_gapo = (opt->max_gapo < 0) ? (int)(0.99 + opt->max_gapo_frac * opt->seed2_length) : opt->max_gapo; 
  d->max_gape = (opt->max_gape < 0) ? (int)(0.99 + opt->max_gape_frac * opt->seed2_length) : opt->max_gape; 

  (*data) = (void*)d;

  return 0;
}

tmap_map_sams_t*
tmap_map1_thread_map_core(void **data, tmap_seq_t *seqs[4], int32_t seq_len,
                          tmap_index_t *index, tmap_map_opt_t *opt)
{
  tmap_map1_thread_data_t *d = (tmap_map1_thread_data_t*)(*data);
  int32_t seed2_len = 0;
  tmap_map_opt_t opt_local = (*opt); // copy over values
  tmap_map_sams_t *sams = NULL;
  tmap_string_t *bases = NULL;

  if((0 < opt->min_seq_len && seq_len < opt->min_seq_len)
     || (0 < opt->max_seq_len && opt->max_seq_len < seq_len)) {
      // go to the next loop
      return tmap_map_sams_init(NULL);
  }

  // not enough bases, ignore
  if(0 < opt->seed_length && seq_len < opt->seed_length){
      return tmap_map_sams_init(NULL);
  }

  if(opt->seed2_length < 0 || seq_len < opt->seed2_length) {
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
  
  // get bases for the reversed sequence
  bases = tmap_seq_get_bases(seqs[2]);	

  // primary width, use seed2 length
  if(d->width_length < seed2_len) {
      d->width_length = seed2_len;
      d->width = tmap_realloc(d->width, (1+d->width_length) * sizeof(tmap_bwt_match_width_t), "d->width");
      memset(d->width, 0, (1+d->width_length) * sizeof(tmap_bwt_match_width_t));
  }
  // NB: use the reversed sequence
  tmap_bwt_match_cal_width_reverse(index->bwt, seed2_len, bases->s + (seq_len - seed2_len), d->width);

  // seed width
  if(0 < opt->seed_length) {
      // NB: use the reversed sequence
      tmap_bwt_match_cal_width_reverse(index->bwt, opt->seed_length, bases->s + (seq_len - opt->seed_length), d->seed_width);
  }

  // NB: use the reverse complimented sequence
  sams = tmap_map1_aux_core(seqs[1], index, d->width, (0 < opt_local.seed_length) ? d->seed_width : NULL, &opt_local, d->stack, seed2_len);

  return sams;
}

tmap_map_sams_t*
tmap_map1_thread_map(void **data, tmap_seq_t **seqs, tmap_index_t *index, tmap_rand_t *rand, tmap_map_opt_t *opt)
{
  int32_t seq_len = 0;;
  tmap_map_sams_t *sams = NULL;

  // sequence length
  seq_len = tmap_seq_get_bases_length(seqs[0]);

  // sequence length not in range
  if((0 < opt->min_seq_len && seq_len < opt->min_seq_len)
     || (0 < opt->max_seq_len && opt->max_seq_len < seq_len)) {
      return tmap_map_sams_init(NULL);
  }

  // not enough bases, ignore
  if(0 < opt->seed_length && seq_len < opt->seed_length){
      return tmap_map_sams_init(NULL);
  }

  // core algorithm; use the reverse
  sams = tmap_map1_thread_map_core(data, seqs, seq_len, index, opt);

  return sams;
}

int32_t
tmap_map1_thread_cleanup(void **data, tmap_map_opt_t *opt)
{
  tmap_map1_thread_data_t *d = (tmap_map1_thread_data_t*)(*data);

  tmap_map1_aux_stack_destroy(d->stack);
  free(d->seed_width);
  free(d->width);
  free(d);
  (*data) = NULL;
  return 0;
}

static void 
tmap_map1_core(tmap_map_driver_t *driver)
{
  // add this algorithm
  tmap_map_driver_add(driver,
                      tmap_map1_init, 
                      tmap_map1_thread_init, 
                      tmap_map1_thread_map, 
                      tmap_map1_thread_cleanup,
                      NULL,
                      driver->opt);
  

  // run the driver
  tmap_map_driver_run(driver);
}

int 
tmap_map1_main(int argc, char *argv[])
{
  tmap_map_driver_t *driver = NULL;
  
  // init
  driver = tmap_map_driver_init(TMAP_MAP_ALGO_MAP1, tmap_map1_mapq);
  driver->opt->algo_stage = 1;

  // get options
  if(1 != tmap_map_opt_parse(argc, argv, driver->opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return tmap_map_opt_usage(driver->opt);
  }
  else { 
      // check command line arguments
      tmap_map_opt_check(driver->opt);
  }

  // run map1
  tmap_map1_core(driver);

  // destroy 
  tmap_map_driver_destroy(driver);

  tmap_progress_print2("terminating successfully");

  return 0;
}
