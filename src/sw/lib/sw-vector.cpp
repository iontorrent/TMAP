/*	$Id: sw-vector.c,v 1.15 2009/06/16 23:26:21 rumble Exp $	*/

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>

#include <mmintrin.h>	/* MMX */
#include <xmmintrin.h>	/* SSE */
#include <emmintrin.h>	/* SSE2 */
//#include <pmmintrin.h>/* SSE3 */

#include <sys/time.h>

//#include "../common/util.h"
#include "sw-vector.h"
//#include "../common/time_counter.h"

#define MAX(a, b) (((a) > (b)) ? (a) : (b))

static int	initialised;
static int8_t  *db, *db_ls, *qr;
static int	dblen, qrlen;
static int16_t *nogap, *b_gap;
static int	a_gap_open, a_gap_ext;
static int	b_gap_open, b_gap_ext;
static int	match, mismatch;

/*
 * Calculate the Smith-Waterman score.
 *
 * This is basically an SSE2 version of Wozniak's vectored implementation, but
 * without a score table. Further, we assume a fixed database and query size,
 * so *nogap and *b_gap must be pre-allocated (the malloc overhead for very
 * small scans is _huge_).
 *
 * NOTE THE FOLLOWING:
 *
 *	1) seqA must be padded with 7 bytes at the beginning and end. The first
 *	   element of seqA should be the first pad byte.
 *
 *	2) seqB must be padded with bytes on the end up to mod 8 characters.
 *	   The first element of seqB should be (of course) the first character.
 *
 *	3) seqA and seqB's padding _must_ be different, otherwise our logic will
 *	   consider the padding as matches!
 *
 *      4) These is no _mm_max_epu16 prior to SSE 4! We must use the signed max
 *         function. Unfortunately, this limits our maximum score to 2^15 - 1, or
 *         32767. Since bad things happen if we roll over, our caller must ensure
 *         that this will not happen.
 */
static int
vect_sw_diff_gap(int8_t *seqA, int lena, int8_t *seqB, int lenb,
                 int8_t *ls_seqA, int initbp, bool is_rna)
{
  int i, j, score = 0;
  __m128i v_score, v_zero, v_match, v_mismatch;
  __m128i v_a_gap_ext, v_a_gap_open_ext;
#ifndef v_b_gap_open_ext
  __m128i v_b_gap_ext, v_b_gap_open_ext;
#endif
  __m128i v_a_gap, v_b_gap, v_nogap;
  __m128i v_last_nogap, v_prev_nogap, v_seq_a, v_seq_b;
  __m128i v_tmp;

  /* shut up icc */
  (void)ls_seqA;
  (void)initbp;

#define SET16(a, e7, e6, e5, e4, e3, e2, e1, e0)      \
  _mm_set_epi16((int16_t)a[e7], (int16_t)a[e6], \
                (int16_t)a[e5], (int16_t)a[e4], \
                (int16_t)a[e3], (int16_t)a[e2], \
                (int16_t)a[e1], (int16_t)a[e0])

  v_score		 = _mm_setzero_si128();
  v_zero		 = _mm_setzero_si128();
  v_match		 = SET16((&match), 0, 0, 0, 0, 0, 0, 0, 0);
  v_mismatch	 = SET16((&mismatch), 0, 0, 0, 0, 0, 0, 0, 0);
  v_a_gap_ext	 = SET16((&a_gap_ext), 0, 0, 0, 0, 0, 0, 0, 0);
  v_a_gap_open_ext = SET16((&a_gap_open), 0, 0, 0, 0, 0, 0, 0, 0);
  v_a_gap_open_ext = _mm_add_epi16(v_a_gap_open_ext, v_a_gap_ext);
  v_b_gap_ext	 = SET16((&b_gap_ext), 0, 0, 0, 0, 0, 0, 0, 0);
  v_b_gap_open_ext = SET16((&b_gap_open), 0, 0, 0, 0, 0, 0, 0, 0);
  v_b_gap_open_ext = _mm_add_epi16(v_b_gap_open_ext, v_b_gap_ext);

  for (i = 0; i < lena + 14; i++) {
      nogap[i] = 0;
      b_gap[i] = (int16_t)-b_gap_open;
  }

  for (i = 0; i < (lenb + 7)/8; i++) {
      int k = i * 8;

      v_b_gap = SET16(b_gap, 6, 6, 5, 4, 3, 2, 1, 0);
      v_nogap = SET16(nogap, 6, 6, 5, 4, 3, 2, 1, 0);
      v_seq_a = SET16(seqA, 0, 0, 1, 2, 3, 4, 5, 6);
      v_seq_b = SET16(seqB, k+7, k+6, k+5, k+4, k+3, k+2, k+1, k+0);

      v_a_gap = v_a_gap_ext;
      v_a_gap = _mm_sub_epi16(v_a_gap, v_a_gap_open_ext);

      v_last_nogap = _mm_setzero_si128();
      v_prev_nogap = _mm_setzero_si128();

      for (j = 0; j < (lena + 7); j++) {
          v_b_gap = _mm_slli_si128(v_b_gap, 2);
          v_b_gap = _mm_insert_epi16(v_b_gap, b_gap[j+7], 0);

          v_nogap = _mm_slli_si128(v_nogap, 2);
          v_nogap = _mm_insert_epi16(v_nogap, nogap[j+7], 0);

          v_seq_a = _mm_slli_si128(v_seq_a, 2);
          v_seq_a = _mm_insert_epi16(v_seq_a, seqA[j+7], 0);

          v_tmp = _mm_sub_epi16(v_last_nogap, v_a_gap_open_ext);
          v_a_gap = _mm_sub_epi16(v_a_gap, v_a_gap_ext);
          v_a_gap = _mm_max_epi16(v_a_gap, v_tmp);

          v_tmp = _mm_sub_epi16(v_nogap, v_b_gap_open_ext);
          v_b_gap = _mm_sub_epi16(v_b_gap, v_b_gap_ext);
          v_b_gap = _mm_max_epi16(v_b_gap, v_tmp);

          /* compute the score (v_last_nogap is a tmp variable) */
          v_last_nogap = _mm_cmpeq_epi16(v_seq_a, v_seq_b);
          v_tmp = _mm_and_si128(v_last_nogap, v_match);
          v_last_nogap = _mm_cmpeq_epi16(v_last_nogap, v_zero);
          v_last_nogap = _mm_and_si128(v_last_nogap, v_mismatch);
          v_tmp = _mm_or_si128(v_tmp, v_last_nogap);

          v_last_nogap = _mm_add_epi16(v_prev_nogap, v_tmp);
          v_last_nogap = _mm_max_epi16(v_last_nogap, v_zero);
          v_last_nogap = _mm_max_epi16(v_last_nogap, v_a_gap);
          v_last_nogap = _mm_max_epi16(v_last_nogap, v_b_gap);

          v_prev_nogap = v_nogap;
          v_nogap = v_last_nogap;

          b_gap[j] = (int16_t)_mm_extract_epi16(v_b_gap, 7);
          nogap[j] = (int16_t)_mm_extract_epi16(v_nogap, 7);

          v_score = _mm_max_epi16(v_score, v_last_nogap);
      }
  }

  /*
   * Ugh. Old gcc can't loop and using _mm_store to an int16_t array
   * breaks strict-aliasing rules.
   */
  assert(score == 0);
  score = MAX(score, _mm_extract_epi16(v_score, 0));
  score = MAX(score, _mm_extract_epi16(v_score, 1));
  score = MAX(score, _mm_extract_epi16(v_score, 2));
  score = MAX(score, _mm_extract_epi16(v_score, 3));
  score = MAX(score, _mm_extract_epi16(v_score, 4));
  score = MAX(score, _mm_extract_epi16(v_score, 5));
  score = MAX(score, _mm_extract_epi16(v_score, 6));
  score = MAX(score, _mm_extract_epi16(v_score, 7));

  return (score);
}

/*
 * Ugly repetition. Yeah, it sucks, but I don't want to macroise the entire
 * thing. In the future, just cut the above, paste and change the function
 * name.
 */

#define v_b_gap_open_ext	v_a_gap_open_ext
#define v_b_gap_ext		v_a_gap_ext

static int
vect_sw_same_gap(int8_t *seqA, int lena, int8_t *seqB, int lenb,
                 int8_t *ls_seqA, int initbp, bool is_rna)
{
  int i, j, score = 0;
  __m128i v_score, v_zero, v_match, v_mismatch;
  __m128i v_a_gap_ext, v_a_gap_open_ext;
#ifndef v_b_gap_open_ext
  __m128i v_b_gap_ext, v_b_gap_open_ext;
#endif
  __m128i v_a_gap, v_b_gap, v_nogap;
  __m128i v_last_nogap, v_prev_nogap, v_seq_a, v_seq_b;
  __m128i v_tmp;

  /* shut up icc */
  (void)ls_seqA;
  (void)initbp;

#define SET16(a, e7, e6, e5, e4, e3, e2, e1, e0)      \
  _mm_set_epi16((int16_t)a[e7], (int16_t)a[e6], \
                (int16_t)a[e5], (int16_t)a[e4], \
                (int16_t)a[e3], (int16_t)a[e2], \
                (int16_t)a[e1], (int16_t)a[e0])

  v_score		 = _mm_setzero_si128();
  v_zero		 = _mm_setzero_si128();
  v_match		 = SET16((&match), 0, 0, 0, 0, 0, 0, 0, 0);
  v_mismatch	 = SET16((&mismatch), 0, 0, 0, 0, 0, 0, 0, 0);
  v_a_gap_ext	 = SET16((&a_gap_ext), 0, 0, 0, 0, 0, 0, 0, 0);
  v_a_gap_open_ext = SET16((&a_gap_open), 0, 0, 0, 0, 0, 0, 0, 0);
  v_a_gap_open_ext = _mm_add_epi16(v_a_gap_open_ext, v_a_gap_ext);
  v_b_gap_ext	 = SET16((&b_gap_ext), 0, 0, 0, 0, 0, 0, 0, 0);
  v_b_gap_open_ext = SET16((&b_gap_open), 0, 0, 0, 0, 0, 0, 0, 0);
  v_b_gap_open_ext = _mm_add_epi16(v_b_gap_open_ext, v_b_gap_ext);

  for (i = 0; i < lena + 14; i++) {
      nogap[i] = 0;
      b_gap[i] = (int16_t)-b_gap_open;
  }

  for (i = 0; i < (lenb + 7)/8; i++) {
      int k = i * 8;

      v_b_gap = SET16(b_gap, 6, 6, 5, 4, 3, 2, 1, 0);
      v_nogap = SET16(nogap, 6, 6, 5, 4, 3, 2, 1, 0);
      v_seq_a = SET16(seqA, 0, 0, 1, 2, 3, 4, 5, 6);
      v_seq_b = SET16(seqB, k+7, k+6, k+5, k+4, k+3, k+2, k+1, k+0);

      v_a_gap = v_a_gap_ext;
      v_a_gap = _mm_sub_epi16(v_a_gap, v_a_gap_open_ext);

      v_last_nogap = _mm_setzero_si128();
      v_prev_nogap = _mm_setzero_si128();

      for (j = 0; j < (lena + 7); j++) {
          v_b_gap = _mm_slli_si128(v_b_gap, 2);
          v_b_gap = _mm_insert_epi16(v_b_gap, b_gap[j+7], 0);

          v_nogap = _mm_slli_si128(v_nogap, 2);
          v_nogap = _mm_insert_epi16(v_nogap, nogap[j+7], 0);

          v_seq_a = _mm_slli_si128(v_seq_a, 2);
          v_seq_a = _mm_insert_epi16(v_seq_a, seqA[j+7], 0);

          v_tmp = _mm_sub_epi16(v_last_nogap, v_a_gap_open_ext);
          v_a_gap = _mm_sub_epi16(v_a_gap, v_a_gap_ext);
          v_a_gap = _mm_max_epi16(v_a_gap, v_tmp);

          v_tmp = _mm_sub_epi16(v_nogap, v_b_gap_open_ext);
          v_b_gap = _mm_sub_epi16(v_b_gap, v_b_gap_ext);
          v_b_gap = _mm_max_epi16(v_b_gap, v_tmp);

          /* compute the score (v_last_nogap is a tmp variable) */
          v_last_nogap = _mm_cmpeq_epi16(v_seq_a, v_seq_b);
          v_tmp = _mm_and_si128(v_last_nogap, v_match);
          v_last_nogap = _mm_cmpeq_epi16(v_last_nogap, v_zero);
          v_last_nogap = _mm_and_si128(v_last_nogap, v_mismatch);
          v_tmp = _mm_or_si128(v_tmp, v_last_nogap);

          v_last_nogap = _mm_add_epi16(v_prev_nogap, v_tmp);
          v_last_nogap = _mm_max_epi16(v_last_nogap, v_zero);
          v_last_nogap = _mm_max_epi16(v_last_nogap, v_a_gap);
          v_last_nogap = _mm_max_epi16(v_last_nogap, v_b_gap);

          v_prev_nogap = v_nogap;
          v_nogap = v_last_nogap;

          b_gap[j] = (int16_t)_mm_extract_epi16(v_b_gap, 7);
          nogap[j] = (int16_t)_mm_extract_epi16(v_nogap, 7);

          v_score = _mm_max_epi16(v_score, v_last_nogap);
      }
  }

  /*
   * Ugh. Old gcc can't loop and using _mm_store to an int16_t array
   * breaks strict-aliasing rules.
   */
  assert(score == 0);
  score = MAX(score, _mm_extract_epi16(v_score, 0));
  score = MAX(score, _mm_extract_epi16(v_score, 1));
  score = MAX(score, _mm_extract_epi16(v_score, 2));
  score = MAX(score, _mm_extract_epi16(v_score, 3));
  score = MAX(score, _mm_extract_epi16(v_score, 4));
  score = MAX(score, _mm_extract_epi16(v_score, 5));
  score = MAX(score, _mm_extract_epi16(v_score, 6));
  score = MAX(score, _mm_extract_epi16(v_score, 7));

  return (score);
}

int
sw_vector(uint8_t *target, int32_t target_len,
          uint8_t *query, int32_t query_len)
{
  int i, score;

  if (!initialised)
    abort();

  assert(target_len > 0 && target_len <= dblen);
  assert(query_len > 0 && query_len <= qrlen);

  memset(db, -1, (dblen + 14) * sizeof(db[0]));
  memset(qr, -2, (qrlen + 14) * sizeof(qr[0]));

  for (i = 0; i < target_len; i++)
    db[i+7] = (int8_t)target[i];

  for (i = 0; i < query_len; i++)
    qr[i+7] = (int8_t)query[i];

  if (a_gap_open == b_gap_open && a_gap_ext == b_gap_ext) {
      score = vect_sw_same_gap(&db[0], target_len, &qr[7], query_len,
                               &db_ls[0], 0, false);
  } else { 
      score = vect_sw_diff_gap(&db[0], target_len, &qr[7], query_len,
                               &db_ls[0], 0, false);
  }

  return (score);
}

int sw_vector_cleanup(void) {
    free(db);
    free(db_ls);
    free(qr);
    free(nogap);
    free(b_gap);
    return 0;
}

int
sw_vector_setup(int _dblen, int _qrlen, int _a_gap_open, int _a_gap_ext,
                int _b_gap_open, int _b_gap_ext, int _match, int _mismatch,
                int _use_colours, bool reset_stats)
{
  if (_match * _qrlen >= 32768) {
      fprintf(stderr, "Error: Match Value is too high/reads are too long. "
              "Please ensure that (Match_Value x your_longest_read_length)"
              " is less than 32768! Try using smaller S-W values.");
      exit(1);
  }

  dblen = _dblen;
  db = (int8_t *)malloc((dblen + 14) * sizeof(db[0]));
  if (db == NULL)
    return (1);

  db_ls = (int8_t *)malloc((dblen + 14) * sizeof(db_ls[0]));
  if (db_ls == NULL)
    return (1);

  qrlen = _qrlen;
  qr = (int8_t *)malloc((qrlen + 14) * sizeof(qr[0]));
  if (qr == NULL)
    return (1);
  nogap = (int16_t *)malloc((dblen + 14) * sizeof(nogap[0]));
  if (nogap == NULL)
    return (1);

  b_gap = (int16_t *)malloc((dblen + 14) * sizeof(b_gap[0]));
  if (b_gap == NULL)
    return (1);

  a_gap_open = -(_a_gap_open);
  a_gap_ext  = -(_a_gap_ext);
  b_gap_open = -(_b_gap_open);
  b_gap_ext  = -(_b_gap_ext);
  match = _match;
  mismatch = _mismatch;

  initialised = 1;

  return (0);
}
