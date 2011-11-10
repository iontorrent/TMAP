/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_error.h"
#include "../../util/tmap_definitions.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../util/tmap_map_util.h"
#include "tmap_map_pairing.h"

static int32_t
tmap_map_pairing_get_strand_diff(tmap_map_sam_t *one, tmap_map_sam_t *two, int32_t strandedness)
{
  if(TMAP_MAP_PAIRING_SAME_STRAND == strandedness) {
      return (one->strand == two->strand) ? 1 : 0;
  }
  else {
      return (one->strand == two->strand) ? 0 : 1;
  }
}

static int32_t
tmap_map_pairing_get_position_diff(tmap_map_sam_t *one, tmap_map_sam_t *two, int32_t positioning, int32_t use_target_end)
{
  uint32_t pos_one, pos_two;
  
  if(1 == use_target_end) {
      pos_one = one->target_end;
      pos_two = two->target_end;
  }
  else {
      pos_one = one->target_start;
      pos_two = two->target_start;
  }

  switch(positioning) {
    case TMAP_MAP_PAIRING_POSITIONING_BA:
      return (pos_one - pos_two);
    case TMAP_MAP_PAIRING_POSITIONING_AB:
    default:
      return (pos_two - pos_one);
  }
}

double
tmap_map_pairing_get_num_std(tmap_map_sam_t *one, tmap_map_sam_t *two, double ins_size_mean, double ins_size_std, double ins_size_std_max_num, 
                            int32_t pen_mm, int32_t strandedness, int32_t positioning, int32_t use_target_end)
{
  int32_t strand_diff, position_diff;
  double num_std;

  //tmap_map_pairing_get_strand_diff(tmap_map_sam_t *one, tmap_map_sam_t *two, int32_t strandedness)
  strand_diff = tmap_map_pairing_get_strand_diff(one, two, strandedness);

  if(0 == strand_diff || one->seqid!= two->seqid) {
      num_std = ins_size_std_max_num;
  }
  else {
      position_diff = tmap_map_pairing_get_position_diff(one, two, positioning, use_target_end);
      num_std = fabs(position_diff - ins_size_mean) / ins_size_std;
      if(ins_size_std_max_num < num_std) num_std = ins_size_std_max_num;
  }
  return num_std;
}

int32_t
tmap_map_pairing_score(tmap_map_sam_t *one, tmap_map_sam_t *two, double ins_size_mean, double ins_size_std, double ins_size_std_max_num, 
                            int32_t pen_mm, int32_t strandedness, int32_t positioning)
{
  double num_std;
  // NB: use the target ends
  num_std = tmap_map_pairing_get_num_std(one, two, ins_size_mean, ins_size_std, ins_size_std_max_num, pen_mm, strandedness, positioning, 1);
  return one->score + two->score + (int)(pen_mm * -1.0 * log10( erfc(M_SQRT1_2 * num_std)) + 0.499);
}

// TODO: have other options
// TODO: fold options into opt
void
tmap_map_pairing_pick_pairs(tmap_map_sams_t *one, tmap_map_sams_t *two, tmap_seq_t *one_seq, tmap_seq_t *two_seq, 
                            double ins_size_mean, double ins_size_std, double ins_size_std_max_num, 
                            int32_t pen_mm, int32_t strandedness, int32_t positioning, tmap_map_opt_t *opt)
{
  int32_t i, j, n_i, n_j;
  int32_t best_score, best_subo_score;
  int32_t n_best, n_best_subo;
  int32_t **scores = NULL; // TODO: make this more efficient
  int32_t mapq;

  n_i = one->n;
  n_j = two->n;
  scores = tmap_malloc(sizeof(int32_t*) * n_i, "scores"); 
  for(i=0;i<n_i;i++) {
      scores[i] = tmap_malloc(sizeof(int32_t) * n_j, "scores[i]");
  }

  // score the pairings
  best_score = best_subo_score = INT32_MIN;
  n_best = n_best_subo = 0;
  for(i=0;i<n_i;i++) {
      for(j=0;j<n_j;j++) {
          // score
          scores[i][j] = tmap_map_pairing_score(&one->sams[i], &two->sams[j], 
                                         ins_size_mean, ins_size_std, ins_size_std_max_num,
                                         pen_mm, strandedness, positioning);
          // update best and next-best
          if(best_score == scores[i][j]) {
              n_best++;
          }
          else if(best_score < scores[i][j]) {
              best_subo_score = best_score;
              n_best_subo = n_best;
              best_score = scores[i][j];
              n_best = 1;
          }
          else if(best_subo_score == scores[i][j]) {
              n_best_subo++;
          }
          else if(best_subo_score < scores[i][j]) {
              best_subo_score = scores[i][j];
              n_best_subo = 1;
          }
      }
  }

  // mapq
  if(best_score < best_subo_score) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }
  if(1 == n_best) {
      mapq = tmap_map_util_mapq_score(tmap_seq_get_bases_length(one_seq) + tmap_seq_get_bases_length(two_seq),
                                      n_best, best_score, n_best_subo, best_subo_score, opt);
      // filter
      // TODO
  }
  else {
      mapq = 0;
      // all should have zero
      // TODO
  }

  // free
  for(i=0;i<n_i;i++) {
      free(scores[i]);
  }
  free(scores);
}
