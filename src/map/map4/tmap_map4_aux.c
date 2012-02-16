#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_vec.h"
#include "../../util/tmap_rand.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_sa.h"
#include "../../index/tmap_index.h"
#include "../../index/tmap_bwt_match.h"
#include "../../index/tmap_bwt_match_hash.h"
#include "../../index/tmap_bwt_smem.h"
#include "../../sw/tmap_sw.h"
#include "../util/tmap_map_util.h"
#include "tmap_map4.h"
#include "tmap_map4_aux.h"

tmap_map4_aux_smem_iter_t *
tmap_map4_aux_smem_iter_init()
{
  tmap_map4_aux_smem_iter_t *iter;
  iter = tmap_calloc(1, sizeof(tmap_map4_aux_smem_iter_t), "iter");
  iter->tmpvec[0] = tmap_calloc(1, sizeof(tmap_bwt_smem_intv_vec_t), "iter->tmpvec[0]");
  iter->tmpvec[1] = tmap_calloc(1, sizeof(tmap_bwt_smem_intv_vec_t), "iter->tmpvec[1]");
  iter->matches   = tmap_calloc(1, sizeof(tmap_bwt_smem_intv_vec_t), "iter->matches");
  return iter;
}

void 
tmap_map4_aux_smem_iter_destroy(tmap_map4_aux_smem_iter_t *iter)
{
  free(iter->tmpvec[0]->a);
  free(iter->tmpvec[1]->a);
  free(iter->tmpvec[0]);
  free(iter->tmpvec[1]);
  free(iter->matches->a);
  free(iter->matches);
  free(iter);
}

static void 
tmap_map4_aux_smem_iter_set_query(tmap_map4_aux_smem_iter_t *iter, int len, const uint8_t *query)
{
  iter->query = query;
  iter->start = 0;
  iter->len = len;
}

static int
tmap_map4_aux_smem_iter_next(tmap_map4_aux_smem_iter_t *iter, tmap_bwt_t *bwt)
{
  iter->tmpvec[0]->n = iter->tmpvec[1]->n = iter->matches->n = 0;
  if (iter->start >= iter->len || iter->start < 0) return -1;
  while (iter->start < iter->len && iter->query[iter->start] > 3) ++iter->start; // skip ambiguous bases
  if (iter->start == iter->len) return -1;
  iter->start = tmap_bwt_smem1(bwt, iter->len, iter->query, iter->start, iter->matches, iter->tmpvec);
  return iter->start;
}

static void
tmap_bwt_smem_intv_copy(tmap_bwt_smem_intv_t *dest, tmap_bwt_smem_intv_t *src)
{
  dest->x[0] = src->x[0];
  dest->x[1] = src->x[1];
  dest->x[2] = src->x[2];
  dest->info = src->info;
}

static tmap_bwt_smem_intv_vec_t *
tmap_bwt_smem_intv_vec_init()
{
  tmap_bwt_smem_intv_vec_t *matches;
  matches = tmap_calloc(1, sizeof(tmap_bwt_smem_intv_vec_t), "matches");
  matches->a = NULL;
  matches->n = matches->m = 0;
  return matches;
}

static void
tmap_bwt_smem_intv_vec_push(tmap_bwt_smem_intv_vec_t *matches, tmap_bwt_smem_intv_t *p)
{
  if (matches->n == matches->m) {
      matches->m = (0 < matches->m) ? (matches->m << 1) : 2;
      matches->a = tmap_realloc(matches->a, sizeof(tmap_bwt_smem_intv_t) * matches->m, "matches->a");
  }
  tmap_bwt_smem_intv_copy(&matches->a[matches->n++], p);
}

static void
tmap_bwt_smem_intv_vec_destroy(tmap_bwt_smem_intv_vec_t *matches)
{
  free(matches->a);
  free(matches);
}

tmap_map_sams_t *
tmap_map4_aux_core(tmap_seq_t *seq,
                   tmap_refseq_t *refseq,
                   tmap_bwt_t *bwt,
                   tmap_sa_t *sa,
                   tmap_bwt_match_hash_t *hash,
                   tmap_map4_aux_smem_iter_t *iter,
                   tmap_rand_t *rand,
                   tmap_map_opt_t *opt)
{
  int32_t i, n;
  int32_t start, by;
  int32_t min_seed_length, max_seed_length;
  tmap_bwt_int_t k;
  tmap_map_sams_t *sams;
  tmap_bwt_smem_intv_vec_t *matches;
  int32_t total = 0;
  int32_t max_repr;
  
  uint8_t *query;
  int32_t query_len;

  sams = tmap_map_sams_init(NULL);

  matches = tmap_bwt_smem_intv_vec_init();

  query = (uint8_t*)tmap_seq_get_bases(seq)->s;
  query_len = tmap_seq_get_bases_length(seq);

  min_seed_length = opt->min_seed_length;
  max_seed_length = opt->max_seed_length;

  if(query_len < min_seed_length) min_seed_length = query_len; 
  if(query_len < max_seed_length) max_seed_length = query_len; 

  start = 0;
  by = (opt->seed_step < 0) ? query_len : opt->seed_step;

  max_repr = opt->max_repr;
  max_repr = (opt->max_iwidth < max_repr) ? opt->max_iwidth : max_repr;
  
  while(start + max_seed_length - 1 < query_len) {
      //fprintf(stderr, "start=%d min_seed_length=%d max_seed_length=%d\n", start, min_seed_length, max_seed_length);
      
      // init iter
      tmap_map4_aux_smem_iter_set_query(iter, max_seed_length, query + start);

      // iterate
      while (0 < tmap_map4_aux_smem_iter_next(iter, bwt)) {
          // go through matches
          for (i = 0; i < iter->matches->n; ++i) {
              tmap_bwt_smem_intv_t *p = &iter->matches->a[i];

              //fprintf(stderr, "EM\t%d\t%d\t%ld\n", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->x[2]);
              // too short
              if ((uint32_t)p->info - (p->info>>32) < min_seed_length) continue;
              // update total
              total++;
              // too many hits?
              if (p->x[2] <= opt->max_iwidth || p->x[2] < max_repr) {
                  // OK
                  tmap_bwt_smem_intv_vec_push(matches, p);
              }
              else if (0 < max_repr) {
                  tmap_bwt_smem_intv_t q;
                  double pr = 0.0;
                  int32_t c;
                  // Keep only representative hits

                  /*
                  p->x[2] = (opt->max_iwidth < max_repr) ? opt->max_iwidth : max_repr;
                  tmap_bwt_smem_intv_vec_push(matches, p);
                  */

                  if (1 == opt->rand_repr) {
                      // choose randomly the representitive hits
                      pr = max_repr / (double)p->x[2];
                      q = *p;
                      for (k = c = 0; k < q.x[2]; ++k) {
                          if (tmap_rand_get(rand) < pr) {
                              // fake
                              p->x[0] = q.x[0] + k;
                              p->x[1] = q.x[1] + k;
                              p->x[2] = 1;
                              // push
                              tmap_bwt_smem_intv_vec_push(matches, p);
                              // update count
                              c++;
                              if(max_repr <= c) break;
                          }
                      }
                      // reset
                      *p = q;
                  }
                  else {
                      // choose uniformly the representitive hits
                      pr = p->x[2] / (double)max_repr;
                      q = *p;
                      for (k = c = 0; k < q.x[2]; ++k, ++c) {
                          if (pr < c) {
                              // fake
                              p->x[0] = q.x[0] + k;
                              p->x[1] = q.x[1] + k;
                              p->x[2] = 1;
                              // push
                              tmap_bwt_smem_intv_vec_push(matches, p);
                              // update count
                              c = 0;
                          }
                      }
                      // reset
                      *p = q;
                  }
              }
          }
      }
      start += by;
  }

  // remove seeds if there were too many repetitive hits
  //fprintf(stderr, "matches->n=%d total=%d\n", (int)matches->n, total);
  if (0 < matches->n && !((matches->n / (double)total) < opt->hit_frac)) {
      // realloc
      for (i = n = 0; i < matches->n; ++i) {
          tmap_bwt_smem_intv_t *p = &matches->a[i];
          n += p->x[2];
      }
      tmap_map_sams_realloc(sams, n);
      // go through the matches
      for (i = n = 0; i < matches->n; ++i) {
          tmap_bwt_smem_intv_t *p = &matches->a[i];
          for (k = 0; k < p->x[2]; ++k) {
              tmap_bwt_int_t pacpos;
              uint32_t seqid, pos;
              uint8_t strand;

              pacpos = p->x[0] + k;
              if(bwt->seq_len < pacpos) pacpos = bwt->seq_len;

              // get the packed position
              pacpos = tmap_sa_pac_pos_hash(sa, bwt, pacpos, hash);
              //fprintf(stderr, "p->x[0]=%llu p->x[1]=%llu p->x[2]=%llu k=%llu bwt->seq_len=%llu pacpos=%llu\n", p->x[0], p->x[1], p->x[2], k, bwt->seq_len, pacpos);

              // convert to reference co-ordinates
              if(0 < tmap_refseq_pac2real(refseq, pacpos, 1, &seqid, &pos, &strand)) {
                  tmap_map_sam_t *s;
                  uint32_t lower, upper, match_length;

                  lower = (uint32_t)p->info;
                  upper = (uint32_t)(p->info >> 32);

                  //fprintf(stderr, "1 seqid:%u pos:%u strand:%d lower=%u upper=%u\n", seqid, pos, strand, lower, upper);
                  // contig boundary
                  if(0 == strand) {
                      match_length = (0 == upper) ? 0 : (upper - 1);
                  }
                  else {
                      //match_length = (query_len < upper) ? 0 : (query_len - upper);
                      match_length = (max_seed_length < upper) ? 0 : (max_seed_length - upper);
                  }
                  if(pos <= match_length) pos = 0;
                  else pos -= match_length;
                  //fprintf(stderr, "2 seqid:%u pos:%u strand:%d lower=%u upper=%u\n", seqid, pos, strand, lower, upper);

                  // save
                  s = &sams->sams[n];

                  // save the hit
                  s->algo_id = TMAP_MAP_ALGO_MAP4;
                  s->algo_stage = opt->algo_stage;
                  s->strand = strand;
                  s->seqid = seqid;
                  s->pos = pos;
                  s->target_len = query_len;
                  if(refseq->annos[seqid].len < s->target_len) {
                      s->target_len = refseq->annos[seqid].len;
                  }
                  s->score_subo = INT32_MIN;

                  // map3 aux data
                  tmap_map_sam_malloc_aux(s);

                  n++;
              }
          }
      }
      // realloc
      tmap_map_sams_realloc(sams, n);
  }

  tmap_bwt_smem_intv_vec_destroy(matches);

  return sams;
}