/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_BWT_CHECK_H_
#define TMAP_BWT_CHECK_H_

// TODO
void 
tmap_bwt_check_core2(tmap_bwt_t *bwt, int32_t length, int32_t print_msg, int32_t print_sa, int32_t warn);

// TODO
void 
tmap_bwt_check_core(const char *fn_fasta, int32_t length, int32_t print_sa, int32_t use_hash);

// TODO
int 
tmap_bwt_check(int argc, char *argv[]);
#endif
