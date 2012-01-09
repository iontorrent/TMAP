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
  pos_one += one->pos;
  pos_two += two->pos;

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

  strand_diff = tmap_map_pairing_get_strand_diff(one, two, strandedness);

  if(0 == strand_diff || one->seqid != two->seqid) {
      num_std = ins_size_std_max_num;
  }
  else {
      position_diff = tmap_map_pairing_get_position_diff(one, two, positioning, use_target_end);
      //if(1 == one->strand) position_diff = 0 - position_diff; // reverse on the opposite strand
      num_std = fabs(position_diff - ins_size_mean) / ins_size_std;
      if(ins_size_std_max_num < num_std) num_std = ins_size_std_max_num;
  }
  return num_std;
}

int32_t
tmap_map_pairing_score(tmap_map_sam_t *one, tmap_map_sam_t *two, double ins_size_mean, double ins_size_std, double ins_size_std_max_num, 
                            int32_t pen_mm, int32_t strandedness, int32_t positioning, uint8_t *proper_pair, double *num_std)
{
  // NB: use the target ends
  (*num_std) = tmap_map_pairing_get_num_std(one, two, ins_size_mean, ins_size_std, ins_size_std_max_num, pen_mm, strandedness, positioning, 1);
  (*proper_pair) = ((*num_std) < ins_size_std_max_num) ? 1 : 0; 
  return one->score + two->score - (int)(pen_mm * -1.0 * log10( erfc(M_SQRT1_2 * (*num_std))) + 0.499);
}

void
tmap_map_pairing_pick_pairs(tmap_map_sams_t *one, tmap_map_sams_t *two, tmap_seq_t *one_seq, tmap_seq_t *two_seq, 
                            tmap_rand_t *rand, tmap_map_opt_t *opt)
{
  int32_t i, j, n_i, n_j;
  int32_t best_score, best_subo_score, best_score_i, best_score_j;
  int32_t n_best, n_best_subo, mapq;
  int32_t **scores = NULL; // TODO: make this more efficient
  uint8_t **proper_pairs = NULL;
  double **num_stds = NULL;

  if(one->n <= 0 || two->n <= 0) return;

  n_i = one->n;
  n_j = two->n;
  scores = tmap_malloc(sizeof(int32_t*) * n_i, "scores"); 
  proper_pairs = tmap_malloc(sizeof(uint8_t*) * n_i, "proper_pairs"); 
  num_stds = tmap_malloc(sizeof(double*) * n_i, "num_stds"); 
  for(i=0;i<n_i;i++) {
      scores[i] = tmap_malloc(sizeof(int32_t) * n_j, "scores[i]");
      proper_pairs[i] = tmap_malloc(sizeof(uint8_t) * n_j, "proper_pairs[i]");
      num_stds[i] = tmap_malloc(sizeof(double) * n_j, "num_stds[i]");
  }

  // reset all the mapping qualities to zero
  for(i=0;i<n_i;i++) {
      one->sams[i].mapq = 0;
      one->sams[i].proper_pair = 0;
  }
  for(j=0;j<n_j;j++) {
      two->sams[j].mapq = 0;
      two->sams[j].proper_pair = 0;
  }

  // score the pairings
  best_score = best_subo_score = INT32_MIN;
  best_score_i = best_score_j = -1;
  n_best = n_best_subo = 0;
  for(i=0;i<n_i;i++) {
      for(j=0;j<n_j;j++) {
          // score
          scores[i][j] = tmap_map_pairing_score(&one->sams[i], &two->sams[j], 
                                         opt->ins_size_mean, opt->ins_size_std, opt->ins_size_std_max_num,
                                         opt->pen_mm, opt->strandedness, opt->positioning, 
                                         &proper_pairs[i][j], &num_stds[i][j]);
          // update best and next-best
          if(best_score == scores[i][j]) {
              n_best++;
          }
          else if(best_score < scores[i][j]) {
              best_score_i = i;
              best_score_j = j;
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
      // update the mapping quality
      one->sams[best_score_i].mapq = two->sams[best_score_j].mapq = mapq;
  }
  else {
      // all should have zero already
  }

  // TODO: adjust suboptimal?
      
  // filter
  if(TMAP_MAP_OPT_ALN_MODE_ALL == opt->aln_output_mode) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }
  else if(TMAP_MAP_OPT_ALN_MODE_ALL == opt->aln_output_mode) {
      // do nothing
  }
  else {
      if(TMAP_MAP_OPT_ALN_MODE_RAND_BEST == opt->aln_output_mode) {
          // re-choose best_score_i and best_score_j
          int32_t r = (int32_t)(tmap_rand_get(rand) * n_best);
          while(0 < r) {
              for(i=0;0<r && i<n_i;i++) {
                  for(j=0;0<r && j<n_j;j++) {
                      if(best_score == scores[i][j]) {
                          best_score_i = i;
                          best_score_j = j;
                          r--;
                      }
                  }
              }
          }
      }
      
      // copy to the front
      if(best_score_i != 0) {
          tmap_map_sam_copy_and_nullify(&one->sams[0], &one->sams[best_score_i]);
      }
      if(best_score_j != 0) {
          tmap_map_sam_copy_and_nullify(&two->sams[0], &two->sams[best_score_j]);
      }

      // HERE
      fprintf(stderr, "best_score=%d best_subo_score=%d\n", best_score, best_subo_score); 
      // re-allocate
      tmap_map_sams_realloc(one, 1);
      tmap_map_sams_realloc(two, 1);
      // paired score
      one->sams[0].pscore = best_score;
      two->sams[0].pscore = best_score;
      // proper pair
      one->sams[0].proper_pair = proper_pairs[best_score_i][best_score_j];
      two->sams[0].proper_pair = proper_pairs[best_score_i][best_score_j];
      // the number of standard deviations from the mean
      one->sams[0].num_stds = num_stds[best_score_i][best_score_j]; 
      two->sams[0].num_stds = num_stds[best_score_i][best_score_j]; 
  }

  // free
  for(i=0;i<n_i;i++) {
      free(scores[i]);
      free(proper_pairs[i]);
      free(num_stds[i]);
  }
  free(scores);
  free(proper_pairs);
  free(num_stds);
}
