/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <config.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#include <unistd.h>
#endif
#include <unistd.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sam_print.h"
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
#include "tmap_map3_aux.h"
#include "tmap_map3.h"

int32_t
tmap_map3_get_seed_length(uint64_t ref_len)
{
  int32_t k = 0;
  while(0 < ref_len) {
      ref_len >>= 2; // divide by two
      k++;
  }
  k += 2;
  return k;
}

static inline void
tmap_map3_mapq(tmap_map_sams_t *sams, int32_t score_thr, int32_t score_match, int32_t aln_output_mode)
{
  int32_t i;
  int32_t n_best = 0;
  int32_t best_score, cur_score, best_subo;
  int32_t n_seeds = 0, tot_seeds = 0;
  int32_t mapq;

  // estimate mapping quality TODO: this needs to be refined
  best_score = INT32_MIN;
  best_subo = INT32_MIN;
  n_best = 0;
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      tot_seeds += sams->sams[i].aux.map3_aux->n_seeds;
      if(best_score < cur_score) {
          best_subo = best_score;
          best_score = cur_score;
          n_best = 1;
          n_seeds = sams->sams[i].aux.map3_aux->n_seeds;
      }
      else if(cur_score == best_score) { // qual
          n_best++;
      }
      else {
          if(best_subo < cur_score) {
              best_subo = cur_score;
          }
          cur_score = sams->sams[i].score_subo;
          if(best_subo < cur_score) {
              best_subo = cur_score;
          } 
      }
  }
  if(1 < n_best) {
      mapq = 0;
  }
  else {
      double c = 0.0;
      c = n_seeds / (double)tot_seeds;
      if(best_subo < score_thr) best_subo = score_thr;
      mapq = (int32_t)(c * (best_score - best_subo) * (250.0 / best_score + 0.03 / score_match) + .499);
      if(mapq > 250) mapq = 250;
      if(mapq <= 0) mapq = 1;
  }
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      if(cur_score == best_score) { // qual
          sams->sams[i].mapq = mapq;
      }
      else {
          sams->sams[i].mapq = 0;
      }
  }
}

// thread data
typedef struct {
    uint8_t *flow_order[2];
} tmap_map3_thread_data_t;

int32_t
tmap_map3_init(tmap_refseq_t *refseq, tmap_map_opt_t *opt)
{
  // adjust opt for opt->score_match
  opt->score_thr *= opt->score_match;

  // set the seed length
  if(-1 == opt->seed_length) {
      opt->seed_length = tmap_map3_get_seed_length(refseq->len);
      tmap_progress_print("setting the seed length to %d", opt->seed_length);
  }


  return 0;
}

int32_t
tmap_map3_thread_init(void **data, tmap_map_opt_t *opt)
{
  tmap_map3_thread_data_t *d = NULL;

  d = tmap_calloc(1, sizeof(tmap_map3_thread_data_t), "d");
  d->flow_order[0] = d->flow_order[1] = NULL;

  (*data) = (void*)d;

  return 0;
}

tmap_map_sams_t*
tmap_map3_thread_map(void **data, tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2], tmap_map_opt_t *opt)
{
  int32_t i;
  tmap_seq_t *seqs[2] = {NULL, NULL};
  tmap_map_sams_t *sams = NULL;
  tmap_map3_thread_data_t *d = (tmap_map3_thread_data_t*)(*data);

  // set up the flow order
  if(0 < opt->hp_diff && NULL == d->flow_order[0]) {
      if(TMAP_SEQ_TYPE_SFF == seq->type) {
          d->flow_order[0] = tmap_malloc(sizeof(uint8_t) * 4, "d->flow_order[0]");
          d->flow_order[1] = tmap_malloc(sizeof(uint8_t) * 4, "d->flow_order[1]");
          for(i=0;i<4;i++) {
              d->flow_order[0][i] = tmap_nt_char_to_int[(int)seq->data.sff->gheader->flow->s[i]];
              d->flow_order[1][i] = 3 - tmap_nt_char_to_int[(int)seq->data.sff->gheader->flow->s[i]];
          }
      }
      else {
          tmap_error("bug encountered", Exit, OutOfRange);
      }
  }

  if((0 < opt->min_seq_len && tmap_seq_get_bases(seq)->l < opt->min_seq_len)
     || (0 < opt->max_seq_len && opt->max_seq_len < tmap_seq_get_bases(seq)->l)) {
      return tmap_map_sams_init();
  }

  // clone the sequence 
  seqs[0] = tmap_seq_clone(seq);
  seqs[1] = tmap_seq_clone(seq);

  // Adjust for SFF
  tmap_seq_remove_key_sequence(seqs[0]);
  tmap_seq_remove_key_sequence(seqs[1]);

  // Adjust for SFF
  tmap_seq_remove_key_sequence(seqs[0]);
  tmap_seq_remove_key_sequence(seqs[1]);

  // reverse compliment
  tmap_seq_reverse_compliment(seqs[1]);

  // convert to integers
  tmap_seq_to_int(seqs[0]);
  tmap_seq_to_int(seqs[1]);

  // align
  sams = tmap_map3_aux_core(seqs, d->flow_order, refseq, bwt[1], sa[1], opt);

  if(-1 == opt->algo_stage) { // not part of mapall
      // remove duplicates
      tmap_map_util_remove_duplicates(sams, opt->dup_window);

      // mapping quality
      tmap_map3_mapq(sams, opt->score_thr, opt->score_match, opt->aln_output_mode);

      // filter the alignments
      tmap_map_sams_filter(sams, opt->aln_output_mode);

      // re-align the alignments in flow-space
      if(TMAP_SEQ_TYPE_SFF == seq->type) {
          tmap_map_util_fsw(seq->data.sff,
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
tmap_map3_thread_cleanup(void **data, tmap_map_opt_t *opt)
{
  tmap_map3_thread_data_t *d = (tmap_map3_thread_data_t*)(*data);

  free(d->flow_order[0]);
  free(d->flow_order[1]);
  free(d);
  (*data) = NULL;
  return 0;
}

static void
tmap_map3_core(tmap_map_opt_t *opt)
{
  // run the driver
  tmap_map_driver_core(tmap_map3_init,
                       tmap_map3_thread_init,
                       tmap_map3_thread_map,
                       tmap_map3_thread_cleanup,
                       opt);
}

int 
tmap_map3_main(int argc, char *argv[])
{
  tmap_map_opt_t *opt = NULL;

  // init opt
  opt = tmap_map_opt_init(TMAP_MAP_ALGO_MAP3);

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

  // run map3
  tmap_map3_core(opt);

  // destroy opt
  tmap_map_opt_destroy(opt);

  tmap_progress_print2("terminating successfully");

  return 0;
}
