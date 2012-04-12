#include <cstring>
#include <sstream>
#include <stdio.h>
#include <iostream>
#define __STDC_LIMIT_MACROS // Seriously, I want these, no
#include <stdint.h>
#include <emmintrin.h>
#include "AffineSWOptimizationHash.h"

using namespace std;

AffineSWOptimizationHash::AffineSWOptimizationHash()
{
  int i;
  size = 1024;
  qsc = qec = -1;
  hash = new hash_t[size];
  for(i=0;i<size;i++) {
      hash[i].tlen = 0;
      hash[i].hash = -1;
      hash[i].dir = 0;
      hash[i].b.clear();
  }
  query = "";
}

static uint64_t hashDNA(const string &s) {
    int j;
    int len = s.size();
    uint64_t h = 1;
    if (len < 16) {
        for(j=0;j<len;j++) h = h * 1337 + s[j];
        return h >> 8;
        //return h;
    }

    __m128i mhash = _mm_set1_epi16(1);
    __m128i mmul = _mm_set1_epi16(13);
    __m128i mzero = _mm_setzero_si128();
    __m128i m0, m1a, m1b;
    int i = 0;
    while (i + 16 < len) {
        m0 = _mm_loadu_si128((__m128i*)&s[i]);
        m1a = _mm_unpacklo_epi8(m0, mzero);
        m1b = _mm_unpackhi_epi8(m0, mzero);
        mhash = _mm_mullo_epi16(mhash, mmul);
        mhash = _mm_add_epi16(mhash, m1a);
        mhash = _mm_mullo_epi16(mhash, mmul);
        mhash = _mm_add_epi16(mhash, m1b);
        i += 16;
    }
    m0 = _mm_loadu_si128((__m128i*)&s[len - 16]);
    m1a = _mm_unpacklo_epi8(m0, mzero);
    m1b = _mm_unpackhi_epi8(m0, mzero);
    mhash = _mm_mullo_epi16(mhash, mmul);
    mhash = _mm_add_epi16(mhash, m1a);
    mhash = _mm_mullo_epi16(mhash, mmul);
    mhash = _mm_add_epi16(mhash, m1b);

    /*
    uint32_t* p = (uint32_t*)&mhash;
    h = h * 1337 + p[0];
    h = h * 1337 + p[1];
    h = h * 1337 + p[2];
    h = h * 1337 + p[3];
    */
    h = h * 1337 + (uint64_t)_mm_extract_epi16(mhash, 0);
    h = h * 1337 + (uint64_t)_mm_extract_epi16(mhash, 1);
    h = h * 1337 + (uint64_t)_mm_extract_epi16(mhash, 2);
    h = h * 1337 + (uint64_t)_mm_extract_epi16(mhash, 3);
    //return h;
    return h >> 8;
}

bool AffineSWOptimizationHash::process(const string &b, const string &a, int _qsc, int _qec,
                                          int mm, int mi, int o, int e, int dir,
                                          int *opt, int *te, int *qe, int *n_best)
{
  int i;
  if(_qsc != qsc || _qec != qec || 0 != a.compare(query)) { // reset the hash?
      // reset the hash
      for(i=0;i<size;i++) {
          hash[i].tlen = 0;
          hash[i].hash = -1;
          hash[i].dir = 0;
          hash[i].b.clear();
      }
      query = a;
      qsc = _qsc;
      qec = _qec;
  }
  else {
      uint32_t hh = hashDNA(b);
      int htpos = hh % size;
      if (hash[htpos].hash == hh && hash[htpos].tlen == (int)b.size() && dir == hash[htpos].dir && 0 == b.compare(hash[htpos].b)) {
          (*opt) = hash[htpos].opt;
          (*te) = hash[htpos].te;
          (*qe) = hash[htpos].qe;
          (*n_best) = hash[htpos].n_best;
          return true;
      }
  }
  return false;
}

void AffineSWOptimizationHash::add(const string &b, const string &a, int _qsc, int _qec,
                                          int mm, int mi, int o, int e, int dir,
                                          int *opt, int *te, int *qe, int *n_best)
{
  uint32_t hh = hashDNA(b);
  int htpos = hh % size;
  hash[htpos].hash = hh;
  hash[htpos].tlen = b.size();
  hash[htpos].opt = (*opt);
  hash[htpos].te = (*te);
  hash[htpos].qe = (*qe);
  hash[htpos].n_best = (*n_best);
  hash[htpos].dir = dir;
  hash[htpos].b = b;
}

AffineSWOptimizationHash::~AffineSWOptimizationHash()
{
  delete [] hash;
}
