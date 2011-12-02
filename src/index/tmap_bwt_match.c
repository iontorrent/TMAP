/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <unistd.h>
#include "tmap_bwt.h"
#include "tmap_bwt_match.h"
#include "tmap_bwt_match_hash.h"

inline void
tmap_bwt_match_occ(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, uint8_t c, tmap_bwt_match_occ_t *next)
{
  tmap_bwt_match_hash_occ(bwt, prev, c, next, NULL);
}

inline void
tmap_bwt_match_2occ(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, uint8_t c, tmap_bwt_match_occ_t *next)
{
  tmap_bwt_match_hash_2occ(bwt, prev, c, next, NULL);
}

inline void
tmap_bwt_match_occ4(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, tmap_bwt_match_occ_t next[4])
{
  tmap_bwt_match_hash_occ4(bwt, prev, next, NULL);
}

inline void
tmap_bwt_match_2occ4(const tmap_bwt_t *bwt, tmap_bwt_match_occ_t *prev, tmap_bwt_match_occ_t next[4])
{
  tmap_bwt_match_hash_2occ4(bwt, prev, next, NULL);
}

inline void
tmap_bwt_match_cal_width_forward(const tmap_bwt_t *bwt, int len, const char *str, tmap_bwt_match_width_t *width)
{
  tmap_bwt_match_hash_cal_width_forward(bwt, len, str, width, NULL);
}

inline void
tmap_bwt_match_cal_width_reverse(const tmap_bwt_t *bwt, int len, const char *str, tmap_bwt_match_width_t *width)
{
  tmap_bwt_match_hash_cal_width_reverse(bwt, len, str, width, NULL);
}

inline uint32_t
tmap_bwt_match_exact(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa)
{
  return tmap_bwt_match_hash_exact(bwt, len, str, match_sa, NULL);
}

inline uint32_t
tmap_bwt_match_exact_reverse(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa)
{
  return tmap_bwt_match_hash_exact_reverse(bwt, len, str, match_sa, NULL);
}

inline uint32_t
tmap_bwt_match_exact_alt(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa)
{
  return tmap_bwt_match_hash_exact_alt(bwt, len, str, match_sa, NULL);
}

inline uint32_t
tmap_bwt_match_exact_alt_reverse(const tmap_bwt_t *bwt, int len, const uint8_t *str, tmap_bwt_match_occ_t *match_sa)
{
  return tmap_bwt_match_hash_exact_alt_reverse(bwt, len, str, match_sa, NULL);
}

inline uint32_t
tmap_bwt_match_invPsi(const tmap_bwt_t *bwt, uint32_t sa_intv, uint32_t k, uint32_t *s)
{
  return tmap_bwt_match_hash_invPsi(bwt, sa_intv, k, s, NULL);
}
