/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_RAND_H
#define TMAP_RAND_H

#include <stdint.h>

#define TMAP_RAND_NN 312
#define TMAP_RAND_MM 156
#define TMAP_RAND_UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define TMAP_RAND_LM 0x7FFFFFFFULL /* Least significant 31 bits */

#define tmap_rand_drand48(_kr) ((tmap_rand_int(_kr) >> 11) * (1.0/9007199254740992.0))
#define tmap_rand_sample(_kr, _k, _cnt) ((*(_cnt))++ < (_k)? *(_cnt) - 1 : tmap_rand_rand(_kr) % *(_cnt))

/*!
  a thread-safe random number generator
  */
typedef struct {
    int mti;
    uint64_t mt[TMAP_RAND_NN];
} tmap_rand_t;

/*!
  @param  seed  the random seed
  @return       the initialized random number generator
 */
tmap_rand_t*
tmap_rand_init(uint64_t seed);

/*!
  @param  r     the random number generator to reinitialize
  @param  seed  the random seed
 */
void
tmap_rand_reinit(tmap_rand_t *r, uint64_t seed);

/*!
  @param  r  the initialized random number generator
  @return    a random 64-bit integer
 */
uint64_t 
tmap_rand_int(tmap_rand_t *r);

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
