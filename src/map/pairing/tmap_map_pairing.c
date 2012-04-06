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

static inline int32_t
tmap_map_pairing_get_left(tmap_map_sam_t *sam, int32_t len)
{
  return sam->pos + sam->result.target_end - len + 1;
}

static inline int32_t
tmap_map_pairing_get_right(tmap_map_sam_t *sam, int32_t len)
{
  return sam->pos + sam->result.target_end;
}

static int32_t
tmap_map_pairing_get_position_diff(tmap_map_sam_t *one, tmap_map_sam_t *two, int32_t one_len, int32_t two_len,
                                   int32_t strandedness, int32_t positioning)
{
  uint32_t pos_one_left, pos_one_right;
  uint32_t pos_two_left, pos_two_right;
  int32_t diff = 0;

  // NB:does not matter about strand, since target end is always computed on the
  // "+" strand
  pos_one_left = tmap_map_pairing_get_left(one, one_len);
  pos_two_left = tmap_map_pairing_get_left(two, two_len);
  pos_one_right = tmap_map_pairing_get_right(one, one_len);
  pos_two_right = tmap_map_pairing_get_right(two, two_len);

  switch(strandedness) { 
    case TMAP_MAP_PAIRING_SAME_STRAND:
      switch(positioning) {
        case TMAP_MAP_PAIRING_POSITIONING_AB:
          diff = (0 == one->strand) ? (pos_two_right - pos_one_left) : (pos_one_right - pos_two_left);
        case TMAP_MAP_PAIRING_POSITIONING_BA:
        default:
          diff = (0 == one->strand) ? (pos_one_right - pos_two_left) : (pos_two_right - pos_one_left);
      }
    case TMAP_MAP_PAIRING_OPPOSITE_STRAND:
    default:
      diff = (0 == one->strand) ? (pos_two_right - pos_one_left) : (pos_one_right - pos_two_left);
  }

  /*
  fprintf(stderr, "one_len=%d one->strand=%u one->pos=%u one->target_end=%u\n", one_len, one->strand, one->pos, one->target_end);
  fprintf(stderr, "two_len=%d two->strand=%u two->pos=%u two->target_end=%u\n", two_len, two->strand, two->pos, two->target_end);
  fprintf(stderr, "pos_one_left=%u pos_one_right=%u\n", pos_one_left, pos_one_right);
  fprintf(stderr, "pos_two_left=%u pos_two_right=%u\n", pos_two_left, pos_two_right);
  fprintf(stderr, "diff=%d\n", diff);
  */
  
  return diff;
}

double
tmap_map_pairing_get_num_std(tmap_map_sam_t *one, tmap_map_sam_t *two, int32_t one_len, int32_t two_len, 
                             double ins_size_mean, double ins_size_std, double ins_size_std_max_num, 
                            int32_t pen_mm, int32_t strandedness, int32_t positioning)
{
  int32_t strand_diff, position_diff;
  double num_std;

  strand_diff = tmap_map_pairing_get_strand_diff(one, two, strandedness);

  if(0 == strand_diff || one->seqid != two->seqid) {
      num_std = ins_size_std_max_num;
  }
  else {
      position_diff = tmap_map_pairing_get_position_diff(one, two, one_len, two_len, strandedness, positioning);
      num_std = fabs(position_diff - ins_size_mean) / ins_size_std;
      if(ins_size_std_max_num < num_std) num_std = ins_size_std_max_num;
  }
  /*
  fprintf(stderr, "strand_diff=%d strandedness=%d num_std=%lf\n", strand_diff, strandedness, num_std);
  fprintf(stderr, "ins_size_mean=%lf\n", ins_size_mean);
  fprintf(stderr, "ins_size_std=%lf\n", ins_size_std);
  fprintf(stderr, "ins_size_std_max_num=%lf\n", ins_size_std_max_num);
  */

  return num_std;
}

int32_t
tmap_map_pairing_score(tmap_map_sam_t *one, tmap_map_sam_t *two, int32_t one_len, int32_t two_len,
                       double ins_size_mean, double ins_size_std, double ins_size_std_max_num, 
                       int32_t pen_mm, int32_t strandedness, int32_t positioning, uint8_t *proper_pair, double *num_std)
{
  // NB: use the target ends
  (*num_std) = tmap_map_pairing_get_num_std(one, two, one_len, two_len, ins_size_mean, ins_size_std, ins_size_std_max_num, pen_mm, strandedness, positioning);
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

  // no proper pairs
  for(i=0;i<n_i;i++) {
      one->sams[i].proper_pair = 0;
  }
  for(j=0;j<n_j;j++) {
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
                                                tmap_seq_get_bases_length(one_seq), tmap_seq_get_bases_length(two_seq),
                                                opt->ins_size_mean, opt->ins_size_std, opt->ins_size_std_max_num,
                                                opt->pen_mm, opt->strandedness, opt->positioning, 
                                                &proper_pairs[i][j], &num_stds[i][j]);
          /*
          fprintf(stderr, "one->pos=%d two->pos=%d one->score=%d two->score=%d scores[i][j]=%d proper_pairs[i][j]=%d num_stds[i][j]=%lf\n",
                  one->sams[i].pos, two->sams[j].pos, 
                  one->sams[i].score, two->sams[j].score, 
                  scores[i][j], proper_pairs[i][j], num_stds[i][j]);
                  */
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
      if(best_score < opt->score_thr) {
          mapq = 0;
      }
      else if(best_subo_score < opt->score_thr) {
          // adjust the suboptimal
          best_subo_score = opt->score_thr;
      }
      mapq = tmap_map_util_mapq_score(tmap_seq_get_bases_length(one_seq) + tmap_seq_get_bases_length(two_seq),
                                      n_best, best_score, n_best_subo, best_subo_score, opt);
  }
  else {
      mapq = 0;
  }
  
  //fprintf(stderr, "best_score=%d best_subo_score=%d mapq=%d\n", best_score, best_subo_score, mapq);

  // TODO: adjust suboptimal?
      
  // filter
  if(TMAP_MAP_OPT_ALN_MODE_ALL == opt->aln_output_mode || TMAP_MAP_OPT_ALN_MODE_ALL_BEST == opt->aln_output_mode) {
      tmap_error("bug encountered", Exit, OutOfRange);
  }
  else {
      if(TMAP_MAP_OPT_ALN_MODE_RAND_BEST == opt->aln_output_mode && 1 < n_best) {
          // re-choose best_score_i and best_score_j
          int32_t r = 1 + (int32_t)(tmap_rand_get(rand) * n_best);
          if(n_best < r) tmap_error("bug encountered", Exit, OutOfRange);
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
      /*
      fprintf(stderr, "best_score_i=%d\n", best_score_i);
      fprintf(stderr, "best_score_j=%d\n", best_score_j);
      fprintf(stderr, "num_stds[best_score_i][best_score_j]=%lf\n", num_stds[best_score_i][best_score_j]);
      */
      one->sams[0].num_stds = num_stds[best_score_i][best_score_j]; 
      two->sams[0].num_stds = num_stds[best_score_i][best_score_j]; 
      // update the mapping quality
      //fprintf(stderr, "mapq=%d n_best=%d proper_pairs[best_score_i][best_score_j]=%d\n", mapq, n_best, proper_pairs[best_score_i][best_score_j]);
      if(1 == n_best) {
          if(1 == proper_pairs[best_score_i][best_score_j]) { // proper pair
              if(0 < one->sams[0].mapq && 0 < two->sams[0].mapq) { // both were chosen in SE mapq analysis
                  // TODO: how to upweight?
                  /*
                  // sum the two mapqs, with no regard to PE mapq
                  int32_t mapq2 = one->sams[0].mapq + two->sams[0].mapq;
                  mapq = (mapq < mapq2) ? mapq : mapq2; // bound by PE mapq
                  one->sams[0].mapq = two->sams[0].mapq = mapq;
                  */
                  // upweight
                  /*
                  one->sams[0].mapq = (250 < mapq + one->sams[0].mapq) ? 250 : one->sams[0].mapq + mapq; 
                  two->sams[0].mapq = (250 < mapq + two->sams[0].mapq) ? 250 : two->sams[0].mapq + mapq; 
                  */
                  if(one->sams[0].mapq < mapq) one->sams[0].mapq = mapq; 
                  if(two->sams[0].mapq < mapq) two->sams[0].mapq = mapq; 
              }
              else { // one of the two had zero SE mapq
                  if(0 < mapq && 0 == one->sams[0].mapq) one->sams[0].mapq = (mapq + 7 < two->sams[0].mapq) ? (mapq + 7) : two->sams[0].mapq;
                  if(0 < mapq && 0 == two->sams[0].mapq) two->sams[0].mapq = (mapq + 7 < one->sams[0].mapq) ? (mapq + 7) : one->sams[0].mapq;
              }
          }
          else { // discordant pair
              if(0 < one->sams[0].mapq && 0 < two->sams[0].mapq) { // both were chosen in SE mapq analysis
                  // downweight
                  mapq = (mapq < 10) ? 0 : mapq - 10;
                  one->sams[0].mapq = (one->sams[0].mapq < 10) ? 0 : one->sams[0].mapq - 10;
                  two->sams[0].mapq = (two->sams[0].mapq < 10) ? 0 : two->sams[0].mapq - 10;
                  /*
                  if(mapq < one->sams[0].mapq) one->sams[0].mapq = mapq;
                  if(mapq < two->sams[0].mapq) two->sams[0].mapq = mapq;
                  */
              }
              else if(0 == one->sams[0].mapq && 0 == two->sams[0].mapq) { // both were ambiguous in SE mapq analysis
                  mapq = (mapq < 20) ? 0 : mapq - 20; // downweight
                  one->sams[0].mapq = two->sams[0].mapq = mapq;
              }
              else if(0 < one->sams[0].mapq) { // one was chosen in SE mapq analysis
                  two->sams[0].mapq = one->sams[0].mapq;
                  if(mapq < two->sams[0].mapq) two->sams[0].mapq = mapq;
              }
              else { // two was chosen in SE mapq analysis
                  one->sams[0].mapq = two->sams[0].mapq;
                  if(mapq < one->sams[0].mapq) one->sams[0].mapq = mapq;
              }
          }
          if(one->sams[0].mapq < 0) one->sams[0].mapq = 1;
          if(two->sams[0].mapq < 0) two->sams[0].mapq = 1;
      }
      else {
          one->sams[0].mapq = two->sams[0].mapq = 0;
      }
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

// NB: rescues two from one
static tmap_map_sams_t*
tmap_map_pairing_read_rescue_helper(tmap_refseq_t *refseq,
                                    tmap_map_sams_t *one, tmap_map_sams_t *two, 
                                    tmap_seq_t *one_seq[2], tmap_seq_t *two_seq[2], 
                                    double ins_size_mean, double ins_size_std,
                                    int32_t strandedness, int32_t positioning, // positioning should be relateive to one/two
                                    int32_t read_rescue_std_num,
                                    tmap_rand_t *rand, tmap_map_opt_t *opt)
{
  int32_t i, j;
  int32_t best, n_best, best_mapq=-1;
  tmap_map_sams_t *sams = NULL;

  sams = tmap_map_sams_init(NULL);
  if(one->n <= 0) return sams;

  best = INT32_MIN;
  n_best = 0;
  for(i=0;i<one->n;i++) {
      if(best == one->sams[i].score) {
          n_best++;
          if(best_mapq != one->sams[i].mapq) {
              tmap_error("bug encountered", Exit, OutOfRange);
          }
      }
      else if(best < one->sams[i].score) {
          best_mapq = one->sams[i].mapq;
          best = one->sams[i].score;
          n_best = 1;
      }
  }

  if(n_best <= 0 || best_mapq < opt->read_rescue_mapq_thr) return sams;

  tmap_map_sams_realloc(sams, n_best);

  // copy over data
  for(i=j=0;i<one->n && 0<n_best;i++) {
      if(best == one->sams[i].score) {
          int32_t one_len = tmap_seq_get_bases_length(one_seq[0]);
          int32_t two_len = tmap_seq_get_bases_length(two_seq[0]);
          int32_t one_left = tmap_map_pairing_get_left(&one->sams[i], one_len);
          int32_t one_right = tmap_map_pairing_get_right(&one->sams[i], one_len);
          int32_t add, sub;
          
          // found a best
          n_best--;

          // copy from one to sams
          tmap_map_sam_copy(&sams->sams[j], &one->sams[i]);

          // infer #2 position from #1
          if(TMAP_MAP_PAIRING_SAME_STRAND == strandedness) {
              if(TMAP_MAP_PAIRING_POSITIONING_AB == positioning) {
                  if(0 == one->sams[i].strand) {
                      add = one_left + ins_size_mean + 1; sub = two_len;
                  }
                  else {
                      add = one_right; sub = ins_size_mean; // genomic forward left-most
                      //add = one_right + two_len; sub = ins_size_mean + 1 // genomic forward right-most
                  }
              }
              else {
                  if(0 == one->sams[i].strand) {
                      add = one_right; sub = ins_size_mean; 
                  }
                  else {
                      add = one_left + ins_size_mean + 1; sub = two_len; // genomic forward left-most
                      //add = one_left + ins_size_mean; sub = 0; // genomic forward right-most
                  }
              }
          }
          else {
              // NB: positioning does not matter on opposite strands
              if(0 == one->sams[i].strand) {
                  add = one_left + ins_size_mean + 1; sub = two_len; // genomic forward left-most
                  //add = one_left + ins_size_mean; sub = 0; // genomic forward right-most
              }
              else {
                  add = one_right; sub = ins_size_mean; // genomic forward left-most
                  //add = one_left + two_len + 1; sub = ins_size_mean; // genomic forward right-most
              }
              // opposite strand
              sams->sams[j].strand = 1 - one->sams[i].strand;
          }

          // check that we do not go off the contig
          if(add <= sub) { // before the beginning
              continue;
          }
          else if(refseq->annos[one->sams[i].seqid].len < add - sub) { // off the end
              continue;
          }
          sams->sams[j].pos = add - sub;

          // this is important for tmap_map_util_sw_gen_score
          sams->sams[j].target_len = two_len; 

          /*
          fprintf(stderr, "RR seqid:%u pos:%u strand:%d\n",
                  sams->sams[j].seqid, sams->sams[j].pos, sams->sams[j].strand);
                  */

          // TODO: annotate as read rescued

          j++;
      }
  }

  if(j < sams->n) {
      tmap_map_sams_realloc(sams, n_best);
      if(0 == sams->n) return sams;
  }

  // generate scores
  // NB: no banding, no seed fraction, (TODO: others?)...
  tmap_map_opt_t opt_local = (*opt);
  opt_local.max_seed_band = 0;
  opt_local.stage_seed_freqc = 0.0;
  opt_local.bw += ins_size_std * read_rescue_std_num;
  sams = tmap_map_util_sw_gen_score(refseq, sams, two_seq, rand, &opt_local);

  return sams;
}

int32_t 
tmap_map_pairing_read_rescue(tmap_refseq_t *refseq, 
                             tmap_map_sams_t *one, tmap_map_sams_t *two, 
                             tmap_seq_t *one_seq[2], tmap_seq_t *two_seq[2], 
                             tmap_rand_t *rand, tmap_map_opt_t *opt)
{
  tmap_map_sams_t *one_rr = NULL, *two_rr = NULL;
  int32_t i, flag = 0;
  
  // Rescue #1 from #2
  one_rr = tmap_map_pairing_read_rescue_helper(refseq, two, one, two_seq, one_seq,
                                               opt->ins_size_mean, opt->ins_size_std,
                                               opt->strandedness, 1-opt->positioning, // NB: update positioning 
                                               opt->read_rescue_std_num,
                                               rand, opt);
  //fprintf(stderr, "RR #1: %d\n", one_rr->n);
  if(0 < one_rr->n) {
      /*
      for(i=0;i<one_rr->n;i++) {
          fprintf(stderr, "Rescued seqid:%u pos:%u score:%d\n",
                  one_rr->sams[i].seqid,
                  one_rr->sams[i].pos,
                  one_rr->sams[i].score);
      }
      */
      i = one->n;
      tmap_map_sams_merge(one, one_rr);
      tmap_map_util_remove_duplicates(one, opt->dup_window, rand);
      if(i < one->n) flag |= 0x1; 
  }
  
  // Rescue #2 from #1
  two_rr = tmap_map_pairing_read_rescue_helper(refseq, one, two, one_seq, two_seq,
                                               opt->ins_size_mean, opt->ins_size_std,
                                               opt->strandedness, opt->positioning, 
                                               opt->read_rescue_std_num,
                                               rand, opt);
  //fprintf(stderr, "RR #2: %d\n", one_rr->n);
  if(0 < two_rr->n) {
      /*
      for(i=0;i<two_rr->n;i++) {
          fprintf(stderr, "Rescued seqid:%u pos:%u score:%d\n",
                  two_rr->sams[i].seqid,
                  two_rr->sams[i].pos,
                  two_rr->sams[i].score);
      }
      */
      i = two->n;
      tmap_map_sams_merge(two, two_rr);
      tmap_map_util_remove_duplicates(two, opt->dup_window, rand);
      if(i < two->n) flag |= 0x2; 
  }

  // free memory
  tmap_map_sams_destroy(one_rr);
  tmap_map_sams_destroy(two_rr);

  return flag;
}
