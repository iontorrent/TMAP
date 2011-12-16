/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
/*

   BWT-Index Construction

   This module constructs BWT and auxiliary data structures.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef TMAP_BWT_GEN_H
#define TMAP_BWT_GEN_H

#include "../util/tmap_definitions.h"

/*! 
  BWT Generation Library
  */

/*! 
  create a bwt FASTA file from a packed FASTA file
  @param  fn_fasta      file name of the FASTA file
  @param  is_large      0 to use the short BWT construction algorith, 1 otherwise (large BWT construction algorithm) 
  @param  occ_interval  the desired occurrence interval
  @param  hash_width    the desired k-mer hash width
  @param  check_hash    1 to validate the hash, 0 otherwise
  */
void 
tmap_bwt_pac2bwt(const char *fn_fasta, uint32_t is_large, int32_t occ_interval, int32_t hash_width, int32_t check_hash);

/*! 
  updates a bwt FASTA file for a new hash width
  @param  fn_fasta      file name of the FASTA file
  @param  hash_width    the desired k-mer hash width
  @param  check_hash    1 to validate the hash, 0 otherwise
  */
void 
tmap_bwt_update_hash(const char *fn_fasta, int32_t hash_width, int32_t check_hash);

/*! 
  @param  T  the input string
  @param  n  the length of the input string
  @return    the primary index if no error occurred, -1 or -2 otherwise
  */
uint32_t 
tmap_bwt_gen_short(uint8_t *T, uint32_t n);

#endif
