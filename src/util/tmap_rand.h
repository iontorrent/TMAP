/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_RAND_H
#define TMAP_RAND_H

// http://en.wikipedia.org/wiki/Multiply-with-carry
#define TMAP_RAND_R_LAG 4096 // must be a power of 2
#define TMAP_RAND_PHI 0x9e3779b9
#define TMAP_RAND_A 18782LL
#define TMAP_RAND_C 362436
#define TMAP_RAND_R 0xfffffffe

/*!
  a thread-safe random number generator
  */
typedef struct {
    uint32_t Q[TMAP_RAND_R_LAG]; /*!< the random seed values */
    uint32_t c; /*!< the carry */
    uint32_t i; /*!< the index */
} tmap_rand_t;

/*!
  @param  seed  the random seed
  @return       the initialized random number generator
 */
tmap_rand_t*
tmap_rand_init(uint32_t seed);

/*!
  @param  r  the initialized random number generator
  @return    a random double between 0 and 1.
 */
double
tmap_rand_get(tmap_rand_t *r);

/*!
  @param  r  the initialized random number generator
*/
void
tmap_rand_destroy(tmap_rand_t *r);

#endif
