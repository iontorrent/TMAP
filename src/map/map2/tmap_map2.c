/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
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
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_sa.h"
#include "../../index/tmap_index.h"
#include "../../io/tmap_seq_io.h"
#include "../../server/tmap_shm.h"
#include "../../sw/tmap_sw.h"
#include "../../sw/tmap_fsw.h"
#include "../util/tmap_map_stats.h"
#include "../util/tmap_map_util.h"
#include "../tmap_map_driver.h"
#include "tmap_map2_mempool.h"
#include "tmap_map2_aux.h"
#include "tmap_map2_core.h"
#include "tmap_map2.h"

// sort by max-score
#define __tmap_map2_sam_sort_score_lt(a, b) ((a).score > (b).score)
TMAP_SORT_INIT(tmap_map2_sam_sort_score, tmap_map_sam_t, __tmap_map2_sam_sort_score_lt)

// thread data
typedef struct {
    tmap_map2_global_mempool_t *pool;
} tmap_map2_thread_data_t;

int32_t
tmap_map2_init(void **data, tmap_refseq_t *refseq, tmap_map_opt_t *opt)
{
  // adjust opt for opt->score_match
  opt->score_thr *= opt->score_match;
  opt->length_coef *= opt->score_match;

  return 0;
}

int32_t
tmap_map2_thread_init(void **data, tmap_map_opt_t *opt)
{
  tmap_map2_thread_data_t *d = NULL;

  d = tmap_calloc(1, sizeof(tmap_map2_thread_data_t), "d");
  d->pool = tmap_map2_global_mempool_init();
  (*data) = (void*)d;

  return 0;
}

tmap_map_sams_t*
tmap_map2_thread_map_core(void **data, tmap_seq_t *seqs[4], int32_t seq_len, tmap_index_t *index, tmap_rand_t *rand, tmap_map_opt_t *opt)
{
  tmap_map_sams_t *sams = NULL;
  tmap_map2_thread_data_t *d = (tmap_map2_thread_data_t*)(*data);

  if((0 < opt->min_seq_len && seq_len < opt->min_seq_len)
     || (0 < opt->max_seq_len && opt->max_seq_len < seq_len)) {
      return tmap_map_sams_init(NULL);
  }

  sams = tmap_map2_aux_core(opt, seqs, index->refseq, index->bwt, index->sa, rand, d->pool);

  return sams;
}

static int32_t
tmap_map2_mapq(tmap_map_sams_t *sams, int32_t seq_len, tmap_map_opt_t *opt)
{
  int32_t i;
  int32_t best_score, best_score_subo, num_best;
  int32_t best_repetitive = 0;

  if(0 == sams->n) {
      return 0;
  }

  // sort by decreasing score
  tmap_sort_introsort(tmap_map2_sam_sort_score, sams->n, sams->sams);

  //Note: assumes that the alignments are sorted by decreasing score
  best_repetitive = (1 == sams->sams[0].aux.map2_aux->flag) ? 1 : 0;
  best_score = sams->sams[0].score;
  num_best = 0;
  for(i=0;i<sams->n;i++) {
      if(0 < i && sams->sams[i-1].score < sams->sams[i].score) { // check assumption
          tmap_error("bug encountered", Exit, OutOfRange);
      }
      if(sams->sams[i].score < sams->sams[0].score) {
          break;
      }
      num_best++;
      if(1 == sams->sams[i].aux.map2_aux->flag) {
          best_repetitive = 1;
      }
      if(0 < sams->sams[i].aux.map2_aux->XI) {
          num_best++; // artificially increase
      }
  }
  // get the best sub-optimal score
  best_score_subo = opt->score_thr;
  for(i=0;i<sams->n;i++) {
      if(best_score_subo < sams->sams[i].score_subo) {
          best_score_subo = sams->sams[i].score_subo;
      }
      if(num_best <= i && best_score_subo < sams->sams[i].score) {
          best_score_subo = sams->sams[i].score;
      }
  }

  if(1 < num_best || 0 < best_repetitive) {
      for(i=0;i<sams->n;i++) {
          sams->sams[i].mapq = 0;
      }
  }
  else if(best_score <= best_score_subo) {
      for(i=0;i<num_best;i++) {
          sams->sams[i].mapq = 1;
      }
      for(;i<sams->n;i++) {
          sams->sams[i].mapq = 0;
      }
  }
  else {
      for(i=0;i<num_best;i++) {
          int32_t qual;
          double c = 1.0;

          if(3 != sams->sams[i].aux.map2_aux->XF) c *= .5;
          if(sams->sams[i].aux.map2_aux->XE < 2) c *= .2;
          qual = (int)(c * (sams->sams[i].score - best_score_subo) * (250.0 / sams->sams[i].score + 0.03 / opt->score_match) + .499);
          if(qual > 250) qual = 250;
          if(sams->sams[i].score == sams->sams[i].score_subo) qual = 0;
          sams->sams[i].mapq = qual;
      }
      for(;i<sams->n;i++) {
          sams->sams[i].mapq = 0;
      }
  }

  return 0;
}

tmap_map_sams_t*
tmap_map2_thread_map(void **data, tmap_seq_t **seqs, tmap_index_t *index, tmap_map_stats_t *stat, tmap_rand_t *rand, tmap_map_opt_t *opt)
{
  tmap_map_sams_t *sams = NULL;
  int32_t seq_len;

  seq_len = tmap_seq_get_bases_length(seqs[0]);
  
  if((0 < opt->min_seq_len && seq_len < opt->min_seq_len)
     || (0 < opt->max_seq_len && opt->max_seq_len < seq_len)) {
      return tmap_map_sams_init(NULL);
  }

  // get the sams
  sams = tmap_map2_thread_map_core(data, seqs, seq_len, index, rand, opt);

  return sams;
}

int32_t
tmap_map2_thread_cleanup(void **data, tmap_map_opt_t *opt)
{
  tmap_map2_thread_data_t *d = (tmap_map2_thread_data_t*)(*data);

  tmap_map2_global_mempool_destroy(d->pool);
  free(d);
  (*data) = NULL;
  return 0;
}

static void
tmap_map2_core(tmap_map_driver_t *driver)
{
  // add this algorithm
  tmap_map_driver_add(driver,
                      tmap_map2_init, 
                      tmap_map2_thread_init, 
                      tmap_map2_thread_map, 
                      tmap_map2_thread_cleanup,
                      NULL,
                      driver->opt);
  

  // run the driver
  tmap_map_driver_run(driver);
}

int 
tmap_map2_main(int argc, char *argv[])
{
  tmap_map_driver_t *driver = NULL;

  // init
  driver = tmap_map_driver_init(TMAP_MAP_ALGO_MAP2, tmap_map2_mapq);
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

  // run map2
  tmap_map2_core(driver);

  // destroy 
  tmap_map_driver_destroy(driver);

  tmap_progress_print2("terminating successfully");

  return 0;
}
