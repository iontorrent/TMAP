/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_INDEX_SIZE_H
#define TMAP_INDEX_SIZE_H

/*! 
  structure to store the command line options for 'tmap indexsize'
  */
typedef struct {
    char *fn_fasta;  /*!< the fasta file name (-f) */
    uint64_t len;  /*!< reference length (-l) */
    int32_t occ_interval;  /*!< the occurrence array interval (-o) */
    int32_t hash_width;  /*!< the occurrence hash width (-w) */
    int32_t sa_interval;  /*!< the suffix array interval (-i) */
    int32_t s;  /*!< use 1024B as 1KB (etc.) (-s) */
} tmap_index_size_opt_t;

/*! 
  main-like function for 'tmap indexsize'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int tmap_index_size(int argc, char *argv[]);

#endif
