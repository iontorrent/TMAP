#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
#include "../util/fmap_alloc.h"
#include "../util/fmap_sort.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwtl.h"
#include "../sw/fmap_sw.h"
#include "fmap_map2_chain.h"
#include "fmap_map2_core.h"
#include "fmap_map2_aux.h"
//#include "utils.h"
//#include "stdaln.h"

#define __left_lt(a, b) ((a).end > (b).end)
FMAP_SORT_INIT(hit, fmap_map2_hit_t, __left_lt)

#define __hitG_lt(a, b) ((a).G > (b).G)
FMAP_SORT_INIT(hitG, fmap_map2_hit_t, __hitG_lt)


int32_t
fmap_map2_aux_resolve_duphits(const fmap_bwt_t *bwt, const fmap_sa_t *sa, fmap_map2_aln_t *b, int32_t IS)
{
  int32_t i, j, n;
  if(b->n == 0) return 0;
  if(NULL != bwt) { // convert to chromosomal coordinates if suitable
      int32_t old_n = b->n;      
      fmap_map2_hit_t *old_hits = b->hits;
      for(i = n = 0; i < b->n; ++i) {
          fmap_map2_hit_t *p = old_hits + i;
          if(p->l - p->k + 1 <= IS) n += p->l - p->k + 1;
          else if(p->G > 0) ++n;
      }
      b->n = b->max = n;
      b->hits = fmap_calloc(b->max, sizeof(fmap_map2_hit_t), "b->hits");
      for(i = j = 0; i < old_n; ++i) {
          fmap_map2_hit_t *p = old_hits + i;
          if(p->l - p->k + 1 <= IS) {
              uint32_t k;
              for(k = p->k; k <= p->l; ++k) {
                  b->hits[j] = *p;
                  b->hits[j].k = fmap_sa_pac_pos(sa, bwt, k);
                  b->hits[j].l = 0;
                  ++j;
              }
          } else if(p->G > 0) {
              b->hits[j] = *p;
              b->hits[j].k = fmap_sa_pac_pos(sa, bwt, p->k);
              b->hits[j].l = 0;
              b->hits[j].flag |= 1;
              ++j;
          }
      }
      free(old_hits);
  }
  fmap_sort_introsort(hitG, b->n, b->hits);
  for(i = 1; i < b->n; ++i) {
      fmap_map2_hit_t *p = b->hits + i;
      if(p->G == 0) break;
      for(j = 0; j < i; ++j) {
          fmap_map2_hit_t *q = b->hits + j;
          int32_t compatible = 1;
          if(q->G == 0) continue;
          if(p->l == 0 && q->l == 0) {
              int32_t qol = (p->end < q->end? p->end : q->end) - (p->beg > q->beg? p->beg : q->beg);
              if(qol < 0) qol = 0;
              if((double)qol / (p->end - p->beg) > FMAP_MAP2_MASK_LEVEL || (double)qol / (q->end - q->beg) > FMAP_MAP2_MASK_LEVEL) {
                  int64_t tol = (int64_t)(p->k + p->len < q->k + q->len? p->k + p->len : q->k + q->len)
                    - (int64_t)(p->k > q->k? p->k : q->k);
                  if((double)tol / p->len > FMAP_MAP2_MASK_LEVEL || (double)tol / q->len > FMAP_MAP2_MASK_LEVEL)
                    compatible = 0;
              }
          }
          if(!compatible) {
              p->G = 0;
              break;
          }
      }
  }
  n = i;
  for(i = j = 0; i < n; ++i) {
      if(b->hits[i].G == 0) continue;
      if(i != j) b->hits[j++] = b->hits[i];
      else ++j;
  }
  b->n = j;
  return b->n;
}

static int32_t 
fmap_map2_aux_resolve_query_overlaps(fmap_map2_aln_t *b, double mask_level)
{
  int32_t i, j, n;
  if(b->n == 0) return 0;
  fmap_sort_introsort(hitG, b->n, b->hits);
  for(i = 1; i < b->n; ++i) {
      fmap_map2_hit_t *p = b->hits + i;
      int32_t all_compatible = 1;
      if(p->G == 0) break;
      for(j = 0; j < i; ++j) {
          fmap_map2_hit_t *q = b->hits + j;
          int64_t tol = 0;
          int32_t qol, compatible = 0;
          double fol;
          if(q->G == 0) continue;
          qol = (p->end < q->end? p->end : q->end) - (p->beg > q->beg? p->beg : q->beg);
          if(qol < 0) qol = 0;
          if(p->l == 0 && q->l == 0) {
              tol = (int64_t)(p->k + p->len < q->k + q->len? p->k + p->len : q->k + q->len)
                - (p->k > q->k? p->k : q->k);
              if(tol < 0) tol = 0;
          }
          fol = (double)qol / (p->end - p->beg < q->end - q->beg? p->end - p->beg : q->end - q->beg);
          if(fol < mask_level || (tol > 0 && qol < p->end - p->beg && qol < q->end - q->beg)) compatible = 1;
          if(!compatible) {
              if(q->G2 < p->G) q->G2 = p->G;
              all_compatible = 0;
          }
      }
      if(!all_compatible) p->G = 0;
  }
  n = i;
  for(i = j = 0; i < n; ++i) {
      if(b->hits[i].G == 0) continue;
      if(i != j) b->hits[j++] = b->hits[i];
      else ++j;
  }
  b->n = j;
  return j;
}

void 
fmap_map2_aln_destroy(fmap_map2_aln_t *a)
{
  int32_t i;
  if(NULL == a) return;
  if(NULL != a->cigar) {
      for(i=0;i<a->n;i++) {
          free(a->cigar[i]);
      }
  }
  free(a->cigar); 
  free(a->n_cigar); 
  free(a->hits);
  free(a);
}

#define __gen_ap(par, opt) do { \
    int32_t i; \
    for(i=0;i<25;i++) { \
        (par).matrix[i] = -(opt)->pen_mm; \
    } \
    for(i=0;i<4;i++) { \
        (par).matrix[i*5+i] = (opt)->score_match; \
    } \
    (par).gap_open = (opt)->pen_gapo; (par).gap_ext = (opt)->pen_gape; \
    (par).gap_end = (opt)->pen_gape; \
    (par).row = 5; (par).band_width = opt->band_width; \
} while(0)

// Note: this is the reverse, not the reverse compliment
#define fmap_map2_rseq_i(_refseq, _i) (fmap_refseq_seq_i(_refseq, _refseq->len-_i-1))

#define fmap_map2_aux_reverse_query(_query, _ql) \
  for(i=0;i<(_ql>>2);i++) { \
      uint8_t tmp = _query[i]; \
      _query[i] = _query[_ql-1-i]; \
      _query[_ql-1-i] = tmp; \
  }

static void 
fmap_map2_aux_extend_left(fmap_map2_opt_t *opt, fmap_map2_aln_t *b, 
                          uint8_t *query,
                          int32_t query_length,
                          fmap_refseq_t *refseq,
                          int32_t is_rev, uint8_t *_mem)
{
  int32_t i, matrix[25];
  uint32_t k;
  uint8_t *target = NULL;
  fmap_sw_param_t par;

  par.matrix = matrix;
  __gen_ap(par, opt);
  // sort according to the descending order of query end
  fmap_sort_introsort(hit, b->n, b->hits);
  target = fmap_calloc(((query_length + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + query_length, sizeof(uint8_t), "target");
  fmap_map2_aux_reverse_query(query, query_length); // reverse the query
  // core loop
  for(i = 0; i < b->n; ++i) {
      fmap_map2_hit_t *p = b->hits + i;
      int32_t lt = ((p->beg + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + query_length;
      int32_t score, j;
      fmap_sw_path_t path;
      p->n_seeds = 1;
      if(p->l || p->k == 0) continue;
      for(j = score = 0; j < i; ++j) {
          fmap_map2_hit_t *q = b->hits + j;
          if(q->beg <= p->beg && q->k <= p->k && q->k + q->len >= p->k + p->len) {
              if(q->n_seeds < (1<<14) - 2) ++q->n_seeds;
              ++score;
          }
      }
      if(score) continue;
      if(lt > p->k) lt = p->k;
      if(is_rev) {
          for(k = p->k - 1, j = 0; k > 0 && j < lt; --k) // FIXME: k=0 not considered!
            target[j++] = fmap_refseq_seq_i(refseq, refseq->len-k-1);
      } else {
          for(k = p->k - 1, j = 0; k > 0 && j < lt; --k) // FIXME: k=0 not considered!
            target[j++] = fmap_refseq_seq_i(refseq, k);
      }
      lt = j;
      score = fmap_sw_extend_core(target, lt, query + query_length - p->beg, p->beg, &par, &path, 0, p->G, _mem);
      if(score > p->G) { // extensible
          p->G = score;
          p->len += path.i;
          p->beg -= path.j;
          p->k -= path.i;
      }
  }
  fmap_map2_aux_reverse_query(query, query_length); // reverse back the query
  free(target);
}

static void 
fmap_map2_aux_extend_right(fmap_map2_opt_t *opt, fmap_map2_aln_t *b, 
                           uint8_t *query,
                           int32_t query_length,
                           fmap_refseq_t *refseq,
                           int32_t is_rev, uint8_t *_mem)
{
  int32_t i, matrix[25];
  uint32_t k;
  uint8_t *target = NULL;
  fmap_sw_param_t par;

  par.matrix = matrix;
  __gen_ap(par, opt);
  target = fmap_calloc(((query_length + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + query_length, sizeof(uint8_t), "target");
  for(i = 0; i < b->n; ++i) {
      fmap_map2_hit_t *p = b->hits + i;
      int32_t lt = ((query_length - p->beg + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + query_length;
      int32_t j, score;
      fmap_sw_path_t path;
      if(p->l) continue;
      if(is_rev) {
          for(k = p->k, j = 0; k < p->k + lt && k < refseq->len; ++k)
            target[j++] = fmap_refseq_seq_i(refseq, refseq->len-k-1);
      } else {
          for(k = p->k, j = 0; k < p->k + lt && k < refseq->len; ++k)
            target[j++] = fmap_refseq_seq_i(refseq, k);
      }
      lt = j;
      score = fmap_sw_extend_core(target, lt, query + p->beg, query_length - p->beg, &par, &path, 0, 1, _mem);
      if(score >= p->G) {
          p->G = score;
          p->len = path.i;
          p->end = path.j + p->beg;
      }
  }
  free(target);
}

/* generate CIGAR array(s) in b->cigar[] */
static void 
fmap_map2_aux_gen_cigar(fmap_map2_opt_t *opt, uint8_t *queries[2], 
                        int32_t query_length, fmap_refseq_t *refseq, fmap_map2_aln_t *b)
{
  uint8_t *target = NULL;
  int32_t i, matrix[25];
  fmap_sw_param_t par;
  fmap_sw_path_t *path;

  par.matrix = matrix;
  __gen_ap(par, opt);
  i = ((query_length + 1) / 2 * opt->score_match + opt->pen_gape) / opt->pen_gape + query_length; // maximum possible target length
  target = fmap_calloc(i, sizeof(uint8_t), "target");
  path = fmap_calloc(i + query_length, sizeof(fmap_sw_path_t), "path");
  // memory clean up for b
  if(b->n < b->max) {
      b->max = b->n;
      b->hits = fmap_realloc(b->hits, b->n * sizeof(fmap_map2_hit_t), "b->hits");
  }
  if(b->cigar) free(b->cigar);
  if(b->n_cigar) free(b->n_cigar);
  b->cigar = (uint32_t**)fmap_calloc(b->max, sizeof(void*), "b->cigar");
  b->n_cigar = (int*)fmap_calloc(b->max, sizeof(int), "b->n_cigar");
  // generate CIGAR
  for(i = 0; i < b->n; ++i) {
      fmap_map2_hit_t *p = b->hits + i;
      uint8_t *query;
      uint32_t k;
      int32_t score, path_len, beg, end;
      if(p->l) continue;
      beg = (p->flag & 0x10)? query_length - p->end : p->beg;
      end = (p->flag & 0x10)? query_length - p->beg : p->end;
      query = queries[(p->flag & 0x10)? 1 : 0] + beg;
      for(k = p->k; k < p->k + p->len; ++k) { // in principle, no out-of-boundary here
        target[k - p->k] = fmap_refseq_seq_i(refseq, k);
      }
      score = fmap_sw_global_core(target, p->len, query, end - beg, &par, path, &path_len);
      b->cigar[i] = fmap_sw_path2cigar(path, path_len, &b->n_cigar[i]);
      if(beg != 0 || end < query_length) { // write soft clipping
          b->cigar[i] = fmap_realloc(b->cigar[i], sizeof(uint32_t) * (b->n_cigar[i] + 2), "b->cigar");
          if(beg != 0) {
              memmove(b->cigar[i] + 1, b->cigar[i], b->n_cigar[i] * 4);
              b->cigar[i][0] = beg<<4 | 4;
              ++b->n_cigar[i];
          }
          if(end < query_length) {
              b->cigar[i][b->n_cigar[i]] = (query_length - end)<<4 | 4;
              ++b->n_cigar[i];
          }
      }
  }
  free(target); free(path);
}

static void 
fmap_map2_aux_merge_hits(fmap_map2_aln_t *b[2], int32_t l, int32_t is_reverse)
{
  int32_t i;
  if(b[0]->n + b[1]->n > b[0]->max) {
      b[0]->max = b[0]->n + b[1]->n;
      b[0]->hits = fmap_realloc(b[0]->hits, b[0]->max * sizeof(fmap_map2_hit_t), "b[0]->hits");
  }
  for(i = 0; i < b[1]->n; ++i) {
      fmap_map2_hit_t *p = b[0]->hits + b[0]->n + i;
      *p = b[1]->hits[i];
      if(is_reverse) {
          int32_t x = p->beg;
          p->beg = l - p->end;
          p->end = l - x;
          p->flag |= 0x10;
      }
  }
  b[0]->n += b[1]->n;
  fmap_map2_aln_destroy(b[1]);
  b[1] = 0;
}

static fmap_map2_aln_t *
fmap_map2_aux_aln(fmap_map2_opt_t *opt, fmap_refseq_t *refseq, 
                  fmap_bwt_t *target_bwt, fmap_sa_t *target_sa,
                  fmap_string_t *seq[2], int32_t is_rev, fmap_map2_global_mempool_t *pool)
{
  fmap_map2_aln_t *b[2], **bb[2];
  int32_t k;
  for(k = 0; k < 2; ++k) {
      fmap_bwtl_t *query = fmap_bwtl_seq2bwtl(seq[k]->l, (uint8_t*)seq[k]->s);
      bb[k] = fmap_map2_core_aln(opt, query, target_bwt, target_sa, pool);
      fmap_bwtl_destroy(query);
  }
  b[0] = bb[0][1]; b[1] = bb[1][1];
  fmap_map2_chain_filter(opt, seq[0]->l, b);
  for(k = 0; k < 2; ++k) {
      fmap_map2_aux_extend_left(opt, bb[k][1], (uint8_t*)seq[k]->s, seq[k]->l, refseq, is_rev, pool->aln_mem);
      fmap_map2_aux_merge_hits(bb[k], seq[k]->l, 0);
      fmap_map2_aux_resolve_duphits(NULL, NULL, bb[k][0], 0);
      fmap_map2_aux_extend_right(opt, bb[k][0], (uint8_t*)seq[k]->s, seq[k]->l, refseq, is_rev, pool->aln_mem);
      b[k] = bb[k][0];
      free(bb[k]);		
  }
  fmap_map2_aux_merge_hits(b, seq[0]->l, 1);
  fmap_map2_aux_resolve_query_overlaps(b[0], opt->mask_level);

  return b[0];
}

/* set ->flag to records the origin of the hit (to forward bwt or reverse bwt) */
static void 
fmap_map2_aux_flag_fr(fmap_map2_aln_t *b[2])
{
  int32_t i, j;
  for(i = 0; i < b[0]->n; ++i) {
      fmap_map2_hit_t *p = b[0]->hits + i;
      p->flag |= 0x10000;
  }
  for(i = 0; i < b[1]->n; ++i) {
      fmap_map2_hit_t *p = b[1]->hits + i;
      p->flag |= 0x20000;
  }
  for(i = 0; i < b[0]->n; ++i) {
      fmap_map2_hit_t *p = b[0]->hits + i;
      for(j = 0; j < b[1]->n; ++j) {
          fmap_map2_hit_t *q = b[1]->hits + i;
          if(q->beg == p->beg && q->end == p->end && q->k == p->k && q->len == p->len && q->G == p->G) {
              q->flag |= 0x30000; p->flag |= 0x30000;
              break;
          }
      }
  }
}

static int32_t 
fmap_map2_aux_fix_cigar(fmap_refseq_t *refseq, fmap_map2_hit_t *p, int32_t n_cigar, uint32_t *cigar)
{
  // FIXME: this routine does not work if the query bridge three reference sequences
  int32_t lq;
  int32_t x, y, i;
  uint32_t coor, seqid, refl;

  fmap_refseq_pac2real(refseq, p->k, p->len, &seqid, &coor);
  refl = refseq->annos[seqid].len;
  x = coor; y = 0;
  // test if the alignment goes beyond the boundary
  for(i = 0; i < n_cigar; ++i) {
      int32_t op = cigar[i]&0xf, ln = cigar[i]>>4;
      if(op == 1 || op == 4 || op == 5) y += ln;
      else if(op == 2) x += ln;
      else x += ln, y += ln;
  }
  lq = y; // length of the query sequence
  if(x > refl) { // then fix it
      int32_t j, nc, mq[2], nlen[2];
      uint32_t *cn, kk = 0;
      nc = mq[0] = mq[1] = nlen[0] = nlen[1] = 0;
      cn = fmap_calloc(n_cigar + 3, sizeof(uint32_t), "cn");
      x = coor; y = 0;
      for(i = j = 0; i < n_cigar; ++i) {
          int32_t op = cigar[i]&0xf, ln = cigar[i]>>4;
          if(op == 4 || op == 5 || op == 1) { // ins or clipping
              y += ln;
              cn[j++] = cigar[i];
          } else if(op == 2) { // del
              if(x + ln >= refl && nc == 0) {
                  cn[j++] = (uint32_t)(lq - y)<<4 | 4;
                  nc = j;
                  cn[j++] = (uint32_t)y<<4 | 4;
                  kk = p->k + (x + ln - refl);
                  nlen[0] = x - coor;
                  nlen[1] = p->len - nlen[0] - ln;
              } else cn[j++] = cigar[i];
              x += ln;
          } else if(op == 0) { // match
              if(x + ln >= refl && nc == 0) {
                  // FIXME: not consider a special case where a split right between M and I
                  cn[j++] = (uint32_t)(refl - x)<<4 | 0; // write M
                  cn[j++] = (uint32_t)(lq - y - (refl - x))<<4 | 4; // write S
                  nc = j;
                  mq[0] += refl - x;
                  cn[j++] = (uint32_t)(y + (refl - x))<<4 | 4;
                  if(x + ln - refl) cn[j++] = (uint32_t)(x + ln - refl)<<4 | 0;
                  mq[1] += x + ln - refl;
                  kk = refseq->annos[seqid].offset + refl;
                  nlen[0] = refl - coor;
                  nlen[1] = p->len - nlen[0];
              } else {
                  cn[j++] = cigar[i];
                  mq[nc?1:0] += ln;
              }
              x += ln; y += ln;
          }
      }
      if(mq[0] > mq[1]) { // then take the first alignment
          n_cigar = nc;
          memcpy(cigar, cn, 4 * nc);
          p->len = nlen[0];
      } else {
          p->k = kk; p->len = nlen[1];
          n_cigar = j - nc;
          memcpy(cigar, cn + nc, 4 * (j - nc));
      }
      free(cn);
  }
  return n_cigar;
}

fmap_map2_sam_t *
fmap_map2_sam_init(int32_t n)
{
  fmap_map2_sam_t *sam = NULL;
  sam = fmap_calloc(n, sizeof(fmap_map2_sam_t), "sams");
  sam->entries = fmap_calloc(n, sizeof(fmap_map2_sam_entry_t), "sams->entries");
  sam->num_entries = n;

  return sam;
}

fmap_map2_sam_t *
fmap_map2_sam_realloc(fmap_map2_sam_t *sam, int32_t n)
{
  int32_t i;

  if(n == sam->num_entries) return sam;
  for(i=n;i<sam->num_entries;i++) { // free if shrinking
      free(sam->entries[i].cigar);
  }
  sam->entries = fmap_realloc(sam->entries, n*sizeof(fmap_map2_sam_entry_t), "sam->entries");
  for(i=sam->num_entries;i<n;i++) { // initialize if expanding
      sam->entries[i].cigar = NULL;
      sam->entries[i].n_cigar = 0;
  }
  sam->num_entries = n;

  return sam;
}

void
fmap_map2_sam_destroy(fmap_map2_sam_t *sam)
{
  int32_t i;
  if(NULL == sam) return;
  for(i=0;i<sam->num_entries;i++) {
      free(sam->entries[i].cigar);
  }
  free(sam->entries);
  free(sam);
}

static fmap_map2_sam_t *
fmap_map1_aux_store_hits(fmap_refseq_t *refseq, fmap_map2_opt_t *opt, 
                         fmap_map2_aln_t *aln)
{
  int32_t i, j;
  fmap_map2_sam_t *sam = NULL;

  if(NULL == aln) return NULL;

  sam = fmap_map2_sam_init(aln->n);

  for(i=j=0;i<aln->n;i++) {
      fmap_map2_hit_t *p = aln->hits + i;
      uint32_t seqid = 0, coor = 0;
      int32_t qual;

      if(p->l == 0) {
          aln->n_cigar[i] = fmap_map2_aux_fix_cigar(refseq, p, aln->n_cigar[i], aln->cigar[i]);
          if(fmap_refseq_pac2real(refseq, p->k, p->len, &seqid, &coor) <= 0) {
              continue; // spanning two or more chromosomes
          }
      }

      sam->entries[j].strand = (p->flag & 0x10) ? 1 : 0; // strand
      sam->entries[j].seqid = seqid;
      sam->entries[j].pos = coor;
      if(p->l == 0) {
          // estimate mapping quality
          double c = 1.0;	
          int32_t subo = p->G2 > opt->score_thr ? p->G2 : opt->score_thr;
          if(p->flag>>16 == 1 || p->flag>>16 == 2) c *= .5;
          if(p->n_seeds < 2) c *= .2;
          qual = (int)(c * (p->G - subo) * (250.0 / p->G + 0.03 / opt->score_match) + .499);
          if(qual > 250) qual = 250;
          if(p->flag & 1) qual = 0;
          sam->entries[j].mapq = qual;

          // copy cigar memory
          sam->entries[j].n_cigar = aln->n_cigar[i];
          sam->entries[j].cigar = aln->cigar[i];
          aln->n_cigar[i] = 0;
          aln->cigar[i] = NULL;
      } 
      else {
          sam->entries[j].mapq = 0;
          sam->entries[j].n_cigar = 0;
          sam->entries[j].cigar = NULL;
      }
      sam->entries[j].AS = p->G;
      sam->entries[j].XS = p->G2;
      sam->entries[j].XF = p->flag>>16;
      sam->entries[j].XE = p->n_seeds;
      if(p->l) {
          sam->entries[j].XI = p->l - p->k + 1;
      }
      else {
          sam->entries[j].XI = 0;
      }
      j++;
  }
  if(j != aln->n) {
      sam = fmap_map2_sam_realloc(sam, j);
  }

  return sam;
}

// TODO: return value
fmap_map2_sam_t *
fmap_map2_aux_core(fmap_map2_opt_t *_opt,
                   fmap_seq_t *query,
                   fmap_refseq_t *refseq,
                   fmap_bwt_t *bwt[2],
                   fmap_sa_t *sa[2],
                   fmap_map2_global_mempool_t *pool)
{
  fmap_map2_opt_t opt;
  fmap_string_t *seq[2]={NULL, NULL};
  fmap_string_t *rseq[2]={NULL, NULL};
  fmap_map2_sam_t *sam = NULL;
  fmap_map2_aln_t *b[2]={NULL,NULL};
  fmap_string_t *bases = NULL;
  uint8_t *_seq[2];
  int32_t i, k, l;

  opt = (*_opt);

  // sequence length
  bases = fmap_seq_get_bases(query);
  l = bases->l;

  // set opt->score_thr
  opt.score_thr = _opt->score_thr; // reset opt->score_thr
  opt.score_thr = _opt->score_thr;
  if(opt.score_thr < log(l) * opt.length_coef) opt.score_thr = (int)(log(l) * opt.length_coef + .499);
  if(pool->max_l < l) { // then enlarge working space for fmap_sw_extend_core()
      int32_t tmp = ((l + 1) / 2 * opt.score_match + opt.pen_gape) / opt.pen_gape + l;
      pool->max_l = l;
      pool->aln_mem = fmap_realloc(pool->aln_mem, sizeof(uint8_t) * (tmp + 2) * 24, "pool->aln_mem");
  }

  // set opt->band_width
  opt.band_width = _opt->band_width;
  k = (l * opt.score_match - 2 * opt.pen_gapo) / (2 * opt.pen_gape + opt.score_match);
  i = (l * opt.score_match - opt.score_match - opt.score_thr) / opt.pen_gape;
  if(k > i) k = i;
  if(k < 1) k = 1; // I do not know if k==0 causes troubles
  opt.band_width= _opt->band_width< k ? _opt->band_width: k;

  // set seq[2] and rseq[2]
  seq[0] = fmap_string_init(l);
  seq[1] = fmap_string_init(l);
  rseq[0] = fmap_string_init(l);
  rseq[1] = fmap_string_init(l);

  // convert sequences to 2-bit representation
  for(i=k=0;i<l;i++) {
      uint8_t c = (uint8_t)nt_char_to_int[(int)bases->s[i]];
      if(c >= 4) { c = (int)(drand48() * 4); ++k; } // FIXME: ambiguous bases are not properly handled
      seq[0]->s[i] = c;
      seq[1]->s[l-1-i] = 3 - c; // reverse compliment
      rseq[0]->s[l-1-i] = c; // reverse 
      rseq[1]->s[i] = 3 - c; // compliment
  }
  seq[0]->l = seq[1]->l = rseq[0]->l = rseq[1]->l = l;

  // score threshold
  if(l - k < opt.score_thr) {
      fmap_string_destroy(seq[0]);
      fmap_string_destroy(seq[1]);
      fmap_string_destroy(rseq[0]);
      fmap_string_destroy(rseq[1]);
      return NULL;
  }

  // alignment
  b[0] = fmap_map2_aux_aln(&opt, refseq, bwt[0], sa[0], seq, 0, pool);
  for(k = 0; k < b[0]->n; ++k)
    if(b[0]->hits[k].n_seeds < opt.seeds_rev) break;
  if(k < b[0]->n) {
      b[1] = fmap_map2_aux_aln(&opt, refseq, bwt[1], sa[1], rseq, 1, pool);
      for(i = 0; i < b[1]->n; ++i) {
          fmap_map2_hit_t *p = b[1]->hits + i;
          int x = p->beg;
          p->beg = l - p->end;
          p->end = l - x;
          if(p->l == 0) p->k = refseq->len - (p->k + p->len);
      }
      fmap_map2_aux_flag_fr(b);
      fmap_map2_aux_merge_hits(b, l, 0);
      fmap_map2_aux_resolve_duphits(NULL, NULL, b[0], 0);
      fmap_map2_aux_resolve_query_overlaps(b[0], opt.mask_level);
  } else b[1] = 0;
  // generate CIGAR and print SAM
  _seq[0] = (uint8_t*)seq[0]->s;
  _seq[1] = (uint8_t*)seq[1]->s;
  fmap_map2_aux_gen_cigar(&opt, _seq, l, refseq, b[0]);
  sam = fmap_map1_aux_store_hits(refseq, &opt, b[0]);
  // free
  free(seq[0]);
  fmap_map2_aln_destroy(b[0]);

  return sam;
}
