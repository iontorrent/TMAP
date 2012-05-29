// Coder: logicmachine
// Submission: 20
// URL: http://community.topcoder.com/longcontest/?module=ViewProblemSolution&pm=11786&rd=15078&cr=22887412&subnum=20
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <xmmintrin.h>
#include <emmintrin.h>
#include "Solution7.h"

using namespace std;

void Solution7::process_with_start_clip(
                                        const string &target, const string &query,
                                        int query_end_clip, int match_score, int mismatch_score,
                                        int gap_open, int gap_extension, int dir, ProcessResult &result)
{
  // problem parameters
  int n = target.size(), m = query.size();

  // answer
  //int opt = -1, position = -1, n_best = 0;
  int opt = -1, qe=-1 , te=-1, n_best = 0;

  // u8 -> i16 bridge parameters
  int i16_begin_line = -1, offset_u8 = 0;

  // ==================== process with u8 ====================
  // get aligned memory block
  const int ceil_n_u8 = (n + 15) & ~15;
  const int block_num_u8 = (n + 15) >> 4;
  BlockU8S *blocks_u8 = malign<BlockU8S, ALIGN_SIZE>(m_workarea_u8);
  uint8_t *X_u8 = reinterpret_cast<uint8_t *>(blocks_u8 + block_num_u8);

    {
      // copy target to aligned buffer
      for(int i = 0; i < ceil_n_u8; ++i){
          int src = (i & 15) * block_num_u8 + (i >> 4);
          blocks_u8[i >> 4].target[i & 15] = (src < n ? target[src] : 0);
      }

      // transform parameters
      int neg_mismatch = -mismatch_score;
      int neg_gap_open = -gap_open;
      int neg_gap_ext = -gap_extension;
      int zero = offset_u8 = 1 + neg_gap_open + neg_gap_ext;

      // setup constant vector
      __m128i v_neg_inf = _mm_setzero_si128();
      __m128i v_pos_inf; v_pos_inf = _mm_cmpeq_epi8(v_pos_inf, v_pos_inf);
      __m128i v_zero; mm_set1_epi8(v_zero, zero);
      __m128i v_mi; mm_set1_epi8(v_mi, neg_mismatch);
      __m128i v_ma; mm_set1_epi8(v_ma, match_score);
      __m128i v_gap_open; mm_set1_epi8(v_gap_open, neg_gap_open);
      __m128i v_gap_ext; mm_set1_epi8(v_gap_ext, neg_gap_ext);

      // clear line buffers
      for(int i = 0; i < block_num_u8; ++i){
          _mm_store_si128((__m128i *)(blocks_u8[i].data[0].M), v_zero);
          _mm_store_si128((__m128i *)(blocks_u8[i].data[0].V), v_neg_inf);
          _mm_store_si128((__m128i *)(blocks_u8[i].data[0].H), v_neg_inf);
      }

      // calculate
      for(int i = 0; i < m; ++i){
          const int load = i & 1, store = 1 - load;
          __m128i v_query; mm_set1_epi8(v_query, query[i]);

          const int tail = block_num_u8 - 1;
          __m128i m0 = _mm_load_si128((__m128i *)(blocks_u8[tail].data[load].M));
          m0 = _mm_or_si128(_mm_slli_si128(m0, 1), _mm_srli_si128(v_zero, 15));
          __m128i v0 = _mm_load_si128((__m128i *)(blocks_u8[tail].data[load].V));
          v0 = _mm_slli_si128(v0, 1);
          __m128i h0 = _mm_load_si128((__m128i *)(blocks_u8[tail].data[load].H));
          h0 = _mm_slli_si128(h0, 1);
          __m128i m2 = v_zero, h2 = v_neg_inf;

          __m128i v_row_max = v_neg_inf;
          for(int j = 0; j < block_num_u8; ++j){
              __m128i v_target = _mm_load_si128((__m128i *)(blocks_u8[j].target));
              __m128i match = _mm_cmpeq_epi8(v_query, v_target);
              __m128i score_add = _mm_and_si128(match, v_ma);
              __m128i score_sub = _mm_andnot_si128(match, v_mi);

              __m128i m3 = _mm_max_epu8(_mm_max_epu8(m0, v0), h0);
              m3 = _mm_adds_epu8(m3, score_add);
              m3 = _mm_subs_epu8(m3, score_sub);
              m3 = _mm_max_epu8(m3, v_zero);
              v_row_max = _mm_max_epu8(v_row_max, m3);
              __m128i m1 = _mm_load_si128((__m128i *)(blocks_u8[j].data[load].M));
              _mm_store_si128((__m128i *)(blocks_u8[j].data[store].M), m3);

              __m128i v1 = _mm_load_si128((__m128i *)(blocks_u8[j].data[load].V));
              __m128i m1o = _mm_subs_epu8(m1, v_gap_open);
              __m128i v3 = _mm_subs_epu8(_mm_max_epu8(m1o, v1), v_gap_ext);
              _mm_store_si128((__m128i *)(blocks_u8[j].data[store].V), v3);

              __m128i m2o = _mm_subs_epu8(m2, v_gap_open);
              __m128i h3 = _mm_subs_epu8(_mm_max_epu8(m2o, h2), v_gap_ext);
              __m128i h1 = _mm_load_si128((__m128i *)(blocks_u8[j].data[load].H));
              _mm_store_si128((__m128i *)(blocks_u8[j].data[store].H), h3);

              m0 = m1; m2 = m3; v0 = v1; h0 = h1; h2 = h3;
          }
          if(_mm_movemask_epi8(_mm_cmpeq_epi8(v_row_max, v_pos_inf))){
              i16_begin_line = i;
              break;
          }
            {
              __m128i m2 = _mm_load_si128((__m128i *)(blocks_u8[tail].data[store].M));
              m2 = _mm_or_si128(_mm_slli_si128(m2, 1), _mm_srli_si128(v_zero, 15));
              __m128i h2 = _mm_load_si128((__m128i *)(blocks_u8[tail].data[store].H));
              h2 = _mm_slli_si128(h2, 1);
              int j = 0;
              while(true){
                  __m128i h3_old = _mm_load_si128((__m128i *)(blocks_u8[j].data[store].H));
                  __m128i m2o = _mm_subs_epu8(m2, v_gap_open);
                  __m128i h3_new = _mm_subs_epu8(_mm_max_epu8(m2o, h2), v_gap_ext);
                  __m128i cmpgt = _mm_cmpeq_epi8(
                                                 _mm_setzero_si128(), _mm_subs_epu8(h3_new, h3_old));
                  int cmp = _mm_movemask_epi8(cmpgt) ^ 0xffff;
                  if(!cmp){ break; }
                  _mm_store_si128((__m128i *)(blocks_u8[j].data[store].H), h3_new);
                  if(j == tail){
                      m2 = _mm_load_si128((__m128i *)(blocks_u8[tail].data[store].M));
                      m2 = _mm_or_si128(_mm_slli_si128(m2, 1), _mm_srli_si128(v_zero, 15));
                      h2 = _mm_load_si128((__m128i *)(blocks_u8[tail].data[store].H));
                      h2 = _mm_slli_si128(h2, 1);
                      j = 0;
                  }else{
                      m2 = _mm_load_si128((__m128i *)(blocks_u8[j].data[store].M));
                      h2 = h3_new;
                      ++j;
                  }
              }
            }
          if(query_end_clip){
              int row_max; mm_reduction_max_epu8(row_max, v_row_max);
              if(row_max >= opt){
                  if(row_max > opt){
                      opt = row_max;
                      //position = dir ? 0 : 0x7fffffff;
                      qe = te = 0;
                      n_best = 0;
                  }
                  __m128i v_opt; mm_set1_epi8(v_opt, opt);
                  __m128i v_count = _mm_setzero_si128();
                  /*
                  if(dir){
                      int local_position = 0;
                      for(int j = 0; j < block_num_u8; ++j){
                          __m128i m = _mm_load_si128((__m128i *)(blocks_u8[j].data[store].M));
                          __m128i f = _mm_cmpeq_epi8(m, v_opt);
                          v_count = _mm_add_epi8(v_count, f);
                          const int mask = _mm_movemask_epi8(f);
                          const int k = 31 - __builtin_clz(mask);
                          const int p = k * block_num_u8 + j;
                          if(mask){ local_position = max(p, local_position); }
                      }
                      local_position |= (i << 16);
                      position = max(position, local_position);
                  }else{
                      int local_position = 0x7fffffff;
                      for(int j = block_num_u8 - 1; j >= 0; --j){
                          __m128i m = _mm_load_si128((__m128i *)(blocks_u8[j].data[store].M));
                          __m128i f = _mm_cmpeq_epi8(m, v_opt);
                          v_count = _mm_add_epi8(v_count, f);
                          const int mask = _mm_movemask_epi8(f);
                          const int k = __builtin_ctz(mask);
                          const int p = k * block_num_u8 + j;
                          if(mask){ local_position = min(p, local_position); }
                      }
                      local_position |= (i << 16);
                      position = min(position, local_position);
                  }
                  */
                  int local_position = 0;
                  for(int j = 0; j < block_num_u8; ++j){
                      __m128i m = _mm_load_si128((__m128i *)(blocks_u8[j].data[store].M));
                      __m128i f = _mm_cmpeq_epi8(m, v_opt);
                      v_count = _mm_add_epi8(v_count, f);
                      const int mask = _mm_movemask_epi8(f);
                      const int k = 31 - __builtin_clz(mask);
                      const int p = k * block_num_u8 + j;
                      if(mask){ local_position = max(p, local_position); }
                  }
                  if(dir) {
                      if(qe < i || (qe == i && te < local_position)) qe = i, te = local_position;
                  } else {
                      if(i < qe || (i == qe && te < local_position)) qe = i, te = local_position;
                  }
                  v_count = _mm_sub_epi8(_mm_setzero_si128(), v_count);
                  v_count = _mm_sad_epu8(_mm_setzero_si128(), v_count);
                  n_best += _mm_extract_epi16(v_count, 0) + _mm_extract_epi16(v_count, 4);
              }
          }else if(i == m - 1){
              const int head_mask_shift = (ceil_n_u8 - n) / block_num_u8;
              const int tail_mask_threshold = block_num_u8 - (ceil_n_u8 - n) % block_num_u8;
              __m128i head_mask; head_mask = _mm_cmpeq_epi8(head_mask, head_mask);
              if(head_mask_shift & 8){ head_mask = _mm_srli_si128(head_mask, 8); }
              if(head_mask_shift & 4){ head_mask = _mm_srli_si128(head_mask, 4); }
              if(head_mask_shift & 2){ head_mask = _mm_srli_si128(head_mask, 2); }
              if(head_mask_shift & 1){ head_mask = _mm_srli_si128(head_mask, 1); }
              __m128i tail_mask = _mm_srli_si128(head_mask, 1);
              v_row_max = v_neg_inf;
              for(int j = 0; j < block_num_u8; ++j){
                  __m128i m = _mm_load_si128((__m128i *)(blocks_u8[j].data[store].M));
                  __m128i v = _mm_load_si128((__m128i *)(blocks_u8[j].data[store].V));
                  __m128i mx = _mm_max_epu8(m, v);
                  if(j < tail_mask_threshold){
                      mx = _mm_and_si128(mx, head_mask);
                  }else{
                      mx = _mm_and_si128(mx, tail_mask);
                  }
                  v_row_max = _mm_max_epu8(v_row_max, mx);
                  _mm_store_si128((__m128i *)(X_u8 + (j << 4)), mx);
              }
              mm_reduction_max_epu8(opt, v_row_max);
              __m128i v_opt; mm_set1_epi8(v_opt, opt);
              __m128i v_count = _mm_setzero_si128();
              /*
              if(dir){
                  position = 0;
                  for(int j = 0; j < block_num_u8; ++j){
                      __m128i m = _mm_load_si128((__m128i *)(X_u8 + (j << 4)));
                      __m128i f = _mm_cmpeq_epi8(m, v_opt);
                      v_count = _mm_add_epi8(v_count, f);
                      const int mask = _mm_movemask_epi8(f);
                      const int k = 31 - __builtin_clz(mask);
                      const int p = k * block_num_u8 + j;
                      if(mask){ position = max(p, position); }
                  }
              }else{
                  position = 0x7fffffff;
                  for(int j = block_num_u8 - 1; j >= 0; --j){
                      __m128i m = _mm_load_si128((__m128i *)(X_u8 + (j << 4)));
                      __m128i f = _mm_cmpeq_epi8(m, v_opt);
                      v_count = _mm_add_epi8(v_count, f);
                      const int mask = _mm_movemask_epi8(f);
                      const int k = __builtin_ctz(mask);
                      const int p = k * block_num_u8 + j;
                      if(mask){ position = min(p, position); }
                  }
              }
              */
              te = 0;
              for(int j = 0; j < block_num_u8; ++j){
                  __m128i m = _mm_load_si128((__m128i *)(X_u8 + (j << 4)));
                  __m128i f = _mm_cmpeq_epi8(m, v_opt);
                  v_count = _mm_add_epi8(v_count, f);
                  const int mask = _mm_movemask_epi8(f);
                  const int k = 31 - __builtin_clz(mask);
                  const int p = k * block_num_u8 + j;
                  if(mask){ te = max(p, te); }
              }
              v_count = _mm_sub_epi8(_mm_setzero_si128(), v_count);
              v_count = _mm_sad_epu8(_mm_setzero_si128(), v_count);
              n_best += _mm_extract_epi16(v_count, 0) + _mm_extract_epi16(v_count, 4);
              qe = i;
              //position |= (i << 16);
          }
      }
      opt -= zero;
    }

  // ==================== process with i16 ====================
  if(i16_begin_line >= 0){
      // get aligned memory block
      const int ceil_n = (n + 7) & ~7;
      const int block_num = (n + 7) >> 3;
      BlockI16S *blocks = malign<BlockI16S, ALIGN_SIZE>(m_workarea_i16);
      int16_t *X = reinterpret_cast<int16_t *>(blocks + block_num);

      // copy target to aligned buffer
      for(int i = 0; i < ceil_n; ++i){
          int src = (i & 7) * block_num + (i >> 3);
          blocks[i >> 3].target[i & 7] = (src < n ? target[src] : 0);
      }

      // constant vectors
      __m128i v_neg_inf; v_neg_inf = _mm_cmpeq_epi16(v_neg_inf, v_neg_inf);
      v_neg_inf = _mm_slli_epi16(v_neg_inf, 15);
      __m128i v_ma; mm_set1_epi16(v_ma, match_score - mismatch_score);
      __m128i v_mi; mm_set1_epi16(v_mi, mismatch_score);
      __m128i v_gap_open; mm_set1_epi16(v_gap_open, gap_open);
      __m128i v_gap_ext; mm_set1_epi16(v_gap_ext, gap_extension);

      // copy line buffers
      const int copy_from = i16_begin_line & 1;
      int dst_block = 0, dst_pos = 0;
      for(int i = 0; i < 16 && dst_pos < 8; ++i){
          for(int j = 0; j < block_num_u8; ++j){
              blocks[dst_block].M[dst_pos] =
                static_cast<int>(blocks_u8[j].data[copy_from].M[i]) - offset_u8;
              blocks[dst_block].V[dst_pos] =
                static_cast<int>(blocks_u8[j].data[copy_from].V[i]) - offset_u8;
              blocks[dst_block].H[dst_pos] =
                static_cast<int>(blocks_u8[j].data[copy_from].H[i]) - offset_u8;
              if(++dst_block == block_num){
                  if(++dst_pos == 8){ break; }
                  dst_block = 0;
              }
          }
      }

      // calculate
      for(int i = i16_begin_line; i < m; ++i){
          __m128i v_query; mm_set1_epi16(v_query, query[i]);

          __m128i m0 = _mm_load_si128((__m128i *)(blocks[block_num - 1].M));
          m0 = _mm_slli_si128(m0, 2);
          __m128i v0 = _mm_load_si128((__m128i *)(blocks[block_num - 1].V));
          v0 = _mm_or_si128(_mm_slli_si128(v0, 2), _mm_srli_si128(v_neg_inf, 14));
          __m128i h0 = _mm_load_si128((__m128i *)(blocks[block_num - 1].H));
          h0 = _mm_or_si128(_mm_slli_si128(h0, 2), _mm_srli_si128(v_neg_inf, 14));
          __m128i m2 = _mm_setzero_si128(), h2 = v_neg_inf;

          __m128i v_row_max = v_neg_inf;
          for(int j = 0; j < block_num; ++j){
              __m128i v_target = _mm_load_si128((__m128i *)(blocks[j].target));
              __m128i match = _mm_cmpeq_epi16(v_query, v_target);
              __m128i score = _mm_add_epi16(_mm_and_si128(match, v_ma), v_mi);

              __m128i m3 = _mm_max_epi16(_mm_max_epi16(m0, v0), h0);
              m3 = _mm_max_epi16(_mm_adds_epi16(m3, score), _mm_setzero_si128());
              if(query_end_clip){ v_row_max = _mm_max_epi16(v_row_max, m3); }
              __m128i m1 = _mm_load_si128((__m128i *)(blocks[j].M));
              _mm_store_si128((__m128i *)(blocks[j].M), m3);

              __m128i v1 = _mm_load_si128((__m128i *)(blocks[j].V));
              __m128i m1o = _mm_adds_epi16(m1, v_gap_open);
              __m128i v3 = _mm_adds_epi16(_mm_max_epi16(m1o, v1), v_gap_ext);
              _mm_store_si128((__m128i *)(blocks[j].V), v3);

              __m128i m2o = _mm_adds_epi16(m2, v_gap_open);
              __m128i h3 = _mm_adds_epi16(_mm_max_epi16(m2o, h2), v_gap_ext);
              __m128i h1 = _mm_load_si128((__m128i *)(blocks[j].H));
              _mm_store_si128((__m128i *)(blocks[j].H), h3);

              m0 = m1; m2 = m3; v0 = v1; h0 = h1; h2 = h3;
          }
            {
              __m128i m2 = _mm_load_si128((__m128i *)(blocks[block_num - 1].M));
              m2 = _mm_slli_si128(m2, 2);
              __m128i h2 = _mm_load_si128((__m128i *)(blocks[block_num - 1].H));
              h2 = _mm_or_si128(_mm_slli_si128(h2, 2), _mm_srli_si128(v_neg_inf, 14));
              int j = 0;
              while(true){
                  __m128i h3_old = _mm_load_si128((__m128i *)(blocks[j].H));
                  __m128i m2o = _mm_adds_epi16(m2, v_gap_open);
                  __m128i h3_new = _mm_adds_epi16(_mm_max_epi16(m2o, h2), v_gap_ext);
                  int cmp = _mm_movemask_epi8(_mm_cmpgt_epi16(h3_new, h3_old));
                  if(!cmp){ break; }
                  _mm_store_si128((__m128i *)(blocks[j].H), h3_new);
                  if(j == block_num - 1){
                      m2 = _mm_load_si128((__m128i *)(blocks[j].M));
                      m2 = _mm_slli_si128(m2, 2);
                      h2 = _mm_load_si128((__m128i *)(blocks[j].H));
                      h2 = _mm_or_si128(_mm_slli_si128(h2, 2), _mm_srli_si128(v_neg_inf, 14));
                      j = 0;
                  }else{
                      m2 = _mm_load_si128((__m128i *)(blocks[j].M));
                      h2 = h3_new;
                      ++j;
                  }
              }
            }
          if(query_end_clip){
              int row_max; mm_reduction_max_epi16(row_max, v_row_max);
              if(row_max >= opt){
                  if(row_max > opt){
                      opt = row_max;
                      //position = dir ? 0 : 0x7fffffff;
                      qe = te = 0;
                      n_best = 0;
                  }
                  __m128i v_opt; mm_set1_epi16(v_opt, opt);
                  __m128i v_count = _mm_setzero_si128();
                  /*
                  if(dir){
                      int local_position = 0;
                      for(int j = 0; j < block_num; ++j){
                          __m128i m = _mm_load_si128((__m128i *)(blocks[j].M));
                          __m128i f = _mm_cmpeq_epi16(m, v_opt);
                          v_count = _mm_add_epi8(v_count, f);
                          const int mask = _mm_movemask_epi8(f);
                          const int k = (31 - __builtin_clz(mask)) >> 1;
                          const int p = k * block_num + j;
                          if(mask){ local_position = max(p, local_position); }
                      }
                      local_position |= (i << 16);
                      position = max(position, local_position);
                  }else{
                      int local_position = 0x7fffffff;
                      for(int j = block_num - 1; j >= 0; --j){
                          __m128i m = _mm_load_si128((__m128i *)(blocks[j].M));
                          __m128i f = _mm_cmpeq_epi16(m, v_opt);
                          v_count = _mm_add_epi8(v_count, f);
                          const int mask = _mm_movemask_epi8(f);
                          const int k = __builtin_ctz(mask) >> 1;
                          const int p = k * block_num + j;
                          if(mask){ local_position = min(p, local_position); }
                      }
                      local_position |= (i << 16);
                      position = min(position, local_position);
                  }
                  */
                  int local_position = 0;
                  for(int j = 0; j < block_num; ++j){
                      __m128i m = _mm_load_si128((__m128i *)(blocks[j].M));
                      __m128i f = _mm_cmpeq_epi16(m, v_opt);
                      v_count = _mm_add_epi8(v_count, f);
                      const int mask = _mm_movemask_epi8(f);
                      const int k = (31 - __builtin_clz(mask)) >> 1;
                      const int p = k * block_num + j;
                      if(mask){ local_position = max(p, local_position); }
                  }
                  if(dir) {
                      if(qe < i || (qe == i && te < local_position)) qe = i, te = local_position;
                  } else {
                      if(i < qe || (i == qe && te < local_position)) qe = i, te = local_position;
                  }
                  v_count = _mm_sub_epi8(_mm_setzero_si128(), v_count);
                  v_count = _mm_sad_epu8(_mm_setzero_si128(), v_count);
                  n_best += (_mm_extract_epi16(v_count, 0) + _mm_extract_epi16(v_count, 4)) >> 1;
              }
          }else if(i == m - 1){
              const int head_mask_shift = (ceil_n - n) / block_num;
              const int tail_mask_threshold = block_num - (ceil_n - n) % block_num;
              __m128i head_mask; head_mask = _mm_cmpeq_epi16(head_mask, head_mask);
              if(head_mask_shift & 4){ head_mask = _mm_srli_si128(head_mask, 8); }
              if(head_mask_shift & 2){ head_mask = _mm_srli_si128(head_mask, 4); }
              if(head_mask_shift & 1){ head_mask = _mm_srli_si128(head_mask, 2); }
              __m128i tail_mask = _mm_srli_si128(head_mask, 2);
              v_row_max = v_neg_inf;
              for(int j = 0; j < block_num; ++j){
                  __m128i m = _mm_load_si128((__m128i *)(blocks[j].M));
                  __m128i v = _mm_load_si128((__m128i *)(blocks[j].V));
                  __m128i mx = _mm_max_epi16(m, v);
                  if(j < tail_mask_threshold){
                      mx = _mm_and_si128(mx, head_mask);
                  }else{
                      mx = _mm_and_si128(mx, tail_mask);
                  }
                  v_row_max = _mm_max_epi16(v_row_max, mx);
                  _mm_store_si128((__m128i *)(X + (j << 3)), mx);
              }
              mm_reduction_max_epi16(opt, v_row_max);
              __m128i v_opt; mm_set1_epi16(v_opt, opt);
              __m128i v_count = _mm_setzero_si128();
              /*
              if(dir){
                  position = 0;
                  for(int j = 0; j < block_num; ++j){
                      __m128i m = _mm_load_si128((__m128i *)(X + (j << 3)));
                      __m128i f = _mm_cmpeq_epi16(m, v_opt);
                      v_count = _mm_add_epi8(v_count, f);
                      const int mask = _mm_movemask_epi8(f);
                      const int k = (31 - __builtin_clz(mask)) >> 1;
                      const int p = k * block_num + j;
                      if(mask){ position = max(p, position); }
                  }
              }else{
                  position = 0x7fffffff;
                  for(int j = block_num - 1; j >= 0; --j){
                      __m128i m = _mm_load_si128((__m128i *)(X + (j << 3)));
                      __m128i f = _mm_cmpeq_epi16(m, v_opt);
                      v_count = _mm_add_epi8(v_count, f);
                      const int mask = _mm_movemask_epi8(f);
                      const int k = __builtin_ctz(mask) >> 1;
                      const int p = k * block_num + j;
                      if(mask){ position = min(p, position); }
                  }
              }
              */
              te = 0;
              for(int j = 0; j < block_num; ++j){
                  __m128i m = _mm_load_si128((__m128i *)(X + (j << 3)));
                  __m128i f = _mm_cmpeq_epi16(m, v_opt);
                  v_count = _mm_add_epi8(v_count, f);
                  const int mask = _mm_movemask_epi8(f);
                  const int k = (31 - __builtin_clz(mask)) >> 1;
                  const int p = k * block_num + j;
                  if(mask){ te = max(p, te); }
              }
              v_count = _mm_sub_epi8(_mm_setzero_si128(), v_count);
              v_count = _mm_sad_epu8(_mm_setzero_si128(), v_count);
              n_best += (_mm_extract_epi16(v_count, 0) + _mm_extract_epi16(v_count, 4)) >> 1;
              qe = i;
              //position |= (i << 16);
          }
      }
  }

  // copy answer to result buffer
  result.opt = opt;
  /*
  result.query_end = position >> 16;
  result.target_end = position & 0xffff;
  */
  result.query_end = qe;
  result.target_end = te;
  result.n_best = n_best;
}

void Solution7::process_without_start_clip(
                                           const string &target, const string &query,
                                           int query_end_clip, int match_score, int mismatch_score,
                                           int gap_open, int gap_extension, int dir,
                                           ProcessResult &result)
{
  // problem parameters
  int n = target.size(), m = query.size();
  int block_num = (n + 15) >> 3;
  int x_size = (n + 7) & ~7;
  int simd_query_size = (m + 15) & ~7;
  const char *c_target = target.c_str();
  const char *c_query = query.c_str();

  // answer
  //int opt = -1, position = -1, n_best = 0;
  int opt = -1, qe = -1, te = -1, n_best = 0;

  // get aligned memory block
  uint8_t *aligned_workarea = malign<uint8_t, ALIGN_SIZE>(m_workarea_i16);
  BlockI16 *blocks = reinterpret_cast<BlockI16 *>(aligned_workarea);
  int16_t *X = reinterpret_cast<int16_t *>(blocks + block_num + 1);
  int16_t *aligned_queries[8];
  aligned_queries[0] = X + x_size;
  for(int i = 1; i < 8; ++i){
      aligned_queries[i] = aligned_queries[i - 1] + simd_query_size;
  }
  ++blocks;

  // copy target to block buffer
  for(int i = 0; i < n; ++i){ blocks[i >> 3].target[i & 7] = c_target[i]; }
  for(int i = n; i & 7; ++i){ blocks[i >> 3].target[i & 7] = -1; }

  // copy query to aligned buffer in reverse order
  for(int i = 0; i < 8; ++i){
      int16_t *line = aligned_queries[i];
      for(int j = 0; j < i; ++j){ line[j] = 0; }
      for(int j = 0; j < m; ++j){ line[i + j] = c_query[m - 1 - j]; }
      for(int j = i + m; j & 7; ++j){ line[j] = 0; }
  }

  // setup constant vector
  __m128i v_neg_inf; v_neg_inf = _mm_cmpeq_epi16(v_neg_inf, v_neg_inf);
  v_neg_inf = _mm_slli_epi16(v_neg_inf, 15);
  __m128i v_mi; mm_set1_epi16(v_mi, mismatch_score);
  __m128i v_ma; mm_set1_epi16(v_ma, match_score - mismatch_score);
  __m128i v_gap_open; mm_set1_epi16(v_gap_open, gap_open);
  __m128i v_gap_ext; mm_set1_epi16(v_gap_ext, gap_extension);

  // clear line buffers
  for(int i = -1; i < block_num; ++i){
      _mm_store_si128((__m128i *)(blocks[i].NM), _mm_setzero_si128());
      _mm_store_si128((__m128i *)(blocks[i].M), _mm_setzero_si128());
      _mm_store_si128((__m128i *)(blocks[i].V), v_neg_inf);
      _mm_store_si128((__m128i *)(blocks[i].H), v_neg_inf);
  }

  // calculate
  for(int i = 0; i < m + n - 1; ++i){
      // fix leftmost value
      blocks[-1].H[7] = gap_open + gap_extension * (i + 1);
      blocks[-1].M[7] = NEG_INF;
      // range calculation
      const int begin = max(0, i - m + 1);
      const int simd_begin = begin & ~7;
      const int block_begin = simd_begin >> 3;
      const int end = min(min(i + 1, begin + m), n);
      const int simd_end = (end - 1) & ~7;
      const int block_end = simd_end >> 3;
      // select query alignment
      const int query_begin = m - 1 - i + simd_begin;
      const int16_t *aligned_query =
        aligned_queries[(8 - (query_begin & 7)) & 7] +
        ((query_begin + 7) & ~7) - simd_begin;

      // prefetch data
      __m128i m2 = _mm_srli_si128(_mm_load_si128((__m128i *)(blocks[block_begin - 1].M)), 14);
      __m128i v2 = _mm_srli_si128(_mm_load_si128((__m128i *)(blocks[block_begin - 1].V)), 14);
      __m128i h2 = _mm_srli_si128(_mm_load_si128((__m128i *)(blocks[block_begin - 1].H)), 14);
      __m128i v_row_max = v_neg_inf;

      // update DP table
      int16_t first_block_maxval[8] __attribute__((aligned(16)));
      for(int j = block_begin; j <= block_end; ++j){
          BlockI16 *block = blocks + j;
          /* compare target and query */
          __m128i q  = _mm_load_si128((__m128i *)(aligned_query + (j << 3)));
          __m128i t  = _mm_load_si128((__m128i *)(block->target));
          __m128i f  = _mm_cmpeq_epi16(q, t);
          __m128i b  = _mm_add_epi16(_mm_and_si128(f, v_ma), v_mi);
          /* calculate M */
          __m128i nm = _mm_load_si128((__m128i *)(block->NM));
          __m128i m4 = _mm_adds_epi16(nm, b);
          /* fetch old M */
          __m128i m3 = _mm_load_si128((__m128i *)(block->M));
          /* write new M */
          _mm_store_si128((__m128i *)(block->M), m4);
          /* calculate V */
          __m128i v3 = _mm_load_si128((__m128i *)(block->V));
          __m128i m3o = _mm_adds_epi16(m3, v_gap_open);
          __m128i v4 = _mm_adds_epi16(_mm_max_epi16(m3o, v3), v_gap_ext);
          _mm_store_si128((__m128i *)(block->V), v4);
          /* calculate H */
          __m128i h3 = _mm_load_si128((__m128i *)(block->H));
          m2 = _mm_or_si128(m2, _mm_slli_si128(m3, 2));
          h2 = _mm_or_si128(h2, _mm_slli_si128(h3, 2));
          __m128i m2o = _mm_adds_epi16(m2, v_gap_open);
          __m128i h4 = _mm_adds_epi16(_mm_max_epi16(m2o, h2), v_gap_ext);
          _mm_store_si128((__m128i *)(block->H), h4);
          /* calculate NM */
          v2 = _mm_or_si128(v2, _mm_slli_si128(v3, 2));
          __m128i next_nm = _mm_max_epi16(_mm_max_epi16(h2, v2), m2);
          _mm_store_si128((__m128i *)(block->NM), next_nm);
          /* rotate m2, v2, h2 */
          m2 = _mm_srli_si128(m3, 14);
          v2 = _mm_srli_si128(v3, 14);
          h2 = _mm_srli_si128(h3, 14);
          /* calculate X */
          if(j == block_begin){
              __m128i x4 = _mm_max_epi16(_mm_max_epi16(m4, v4), h4);
              v_row_max = x4;
              _mm_store_si128((__m128i *)first_block_maxval, x4);
          }else if(query_end_clip){
              v_row_max = _mm_max_epi16(v_row_max, m4);
          }
      }

      // fix rightmost value
      BlockI16 *end_block = blocks + (end >> 3);
      end_block->M[end & 7] = 0;
      end_block->V[end & 7] = end_block->H[end & 7] = NEG_INF;

      // update answer
      if(query_end_clip){
          int row_max; mm_reduction_max_epi16(row_max, v_row_max);
          if(row_max >= opt){
              if(row_max > opt){
                  opt = row_max;
                  //position = dir ? 0 : 0x7fffffff;
                  qe = te = 0;
                  n_best = 0;
              }
              mm_set1_epi16(v_row_max, opt);
              __m128i v_best_num = _mm_setzero_si128();
              __m128i temporary_m0;
              if(block_begin == 0){
                  temporary_m0 = _mm_load_si128((__m128i *)(blocks[0].M));
                  __m128i first_block = _mm_load_si128((__m128i *)(first_block_maxval));
                  _mm_store_si128((__m128i *)(blocks[0].M), first_block);
              }
              /*
              if(dir){
                  int local_position = 0;
                  for(int j = block_end; j >= block_begin; --j){
                      __m128i x = _mm_load_si128((__m128i *)(blocks[j].M));
                      __m128i f = _mm_cmpeq_epi16(v_row_max, x);
                      v_best_num = _mm_add_epi8(v_best_num, f);
                      const int mask = _mm_movemask_epi8(f);
                      const int p = (j << 3) + (__builtin_ctz(mask) >> 1);
                      const int pos = ((i - p) << 16) | p;
                      local_position = mask ? pos : local_position;
                  }
                  position = max(position, local_position);
              }else{
                  int local_position = 0;
                  for(int j = block_begin; j <= block_end; ++j){
                      __m128i x = _mm_load_si128((__m128i *)(blocks[j].M));
                      __m128i f = _mm_cmpeq_epi16(v_row_max, x);
                      v_best_num = _mm_add_epi8(v_best_num, f);
                      const int mask = _mm_movemask_epi8(f);
                      const int p = (j << 3) + ((31 - __builtin_clz(mask)) >> 1);
                      const int pos = ((i - p) << 16) | p;
                      local_position = mask ? pos : local_position;
                  }
                  position = min(position, local_position);
              }
              */
              int local_position = 0;
              for(int j = block_end; j >= block_begin; --j){
                  __m128i x = _mm_load_si128((__m128i *)(blocks[j].M));
                  __m128i f = _mm_cmpeq_epi16(v_row_max, x);
                  v_best_num = _mm_add_epi8(v_best_num, f);
                  const int mask = _mm_movemask_epi8(f);
                  const int p = (j << 3) + (__builtin_ctz(mask) >> 1);
                  const int pos = ((i - p) << 16) | p;
                  local_position = mask ? pos : local_position;
              }
              te = max(te, local_position);
              __m128i m = _mm_sub_epi8(_mm_setzero_si128(), v_best_num);
              __m128i sad = _mm_sad_epu8(_mm_setzero_si128(), m);
              n_best += (_mm_extract_epi16(sad, 0) + _mm_extract_epi16(sad, 4)) >> 1;
              if(block_begin == 0){
                  _mm_store_si128((__m128i *)(blocks[0].M), temporary_m0);
              }
          }
      }else{
          X[begin] = first_block_maxval[begin & 7];
      }
  }
  if(!query_end_clip){
      const int block_begin = 0;
      const int block_end = (n + 7) >> 3;
      __m128i v_opt = v_neg_inf;
      for(int i = n; i & 7; ++i){ X[i] = NEG_INF; }
      for(int i = block_begin; i < block_end; ++i){
          __m128i x = _mm_load_si128((__m128i *)(X + (i << 3)));
          v_opt = _mm_max_epi16(v_opt, x);
      }
      mm_reduction_max_epi16(opt, v_opt);
      mm_set1_epi16(v_opt, opt);
      __m128i v_count = _mm_setzero_si128();
      /*
      if(dir){
          for(int i = block_begin; i < block_end; ++i){
              __m128i x = _mm_load_si128((__m128i *)(X + (i << 3)));
              __m128i f = _mm_cmpeq_epi16(x, v_opt);
              v_count = _mm_add_epi8(v_count, f);
              const int mask = _mm_movemask_epi8(f);
              const int k = (31 - __builtin_clz(mask)) >> 1;
              if(mask){ position = (i << 3) + k; }
          }
      }else{
          for(int i = block_end - 1; i >= 0; --i){
              __m128i x = _mm_load_si128((__m128i *)(X + (i << 3)));
              __m128i f = _mm_cmpeq_epi16(x, v_opt);
              v_count = _mm_add_epi8(v_count, f);
              const int mask = _mm_movemask_epi8(f);
              const int k = __builtin_ctz(mask) >> 1;
              if(mask){ position = (i << 3) + k; }
          }
      }
      */
      for(int i = block_begin; i < block_end; ++i){
          __m128i x = _mm_load_si128((__m128i *)(X + (i << 3)));
          __m128i f = _mm_cmpeq_epi16(x, v_opt);
          v_count = _mm_add_epi8(v_count, f);
          const int mask = _mm_movemask_epi8(f);
          const int k = (31 - __builtin_clz(mask)) >> 1;
          if(mask){ te = (i << 3) + k; }
      }
      v_count = _mm_sub_epi8(_mm_setzero_si128(), v_count);
      v_count = _mm_sad_epu8(_mm_setzero_si128(), v_count);
      n_best = (_mm_extract_epi16(v_count, 0) + _mm_extract_epi16(v_count, 4)) >> 1;
      //position |= (m - 1) << 16;
      qe = (m - 1) << 16;
  }

  result.opt = opt;
  /*
  result.query_end = position >> 16;
  result.target_end = position & 0xffff;
  */
  result.query_end = qe;
  result.target_end = te;
  result.n_best = n_best;
}

int Solution7::process(
                          const string& target, const string& query,
                          int query_start_clip, int query_end_clip,
                          int match_score, int mismatch_score,
                          int gap_open, int gap_extension, int dir,
                          int *opt, int *te, int *qe, int *n_best)
{
  ProcessResult result = { 0, 0, 0, 0 };

  if(query_start_clip){
      process_with_start_clip(
                              target, query, query_end_clip,
                              match_score, mismatch_score, gap_open, gap_extension, dir, result);
  }else{
      process_without_start_clip(
                                 target, query, query_end_clip,
                                 match_score, mismatch_score, gap_open, gap_extension, dir, result);
  }

  /*
  ostringstream oss;
  oss << result.opt << " " << result.query_end << " "
    << result.target_end << " " << result.n_best;
  return oss.str();
  */
  (*opt) = result.opt;
  (*te) = result.target_end;
  (*qe) = result.query_end;
  (*n_best) = result.n_best;
  return result.opt;
}

Solution7::Solution7() 
{
  max_qlen = 512;
  max_tlen = 1024;
}
