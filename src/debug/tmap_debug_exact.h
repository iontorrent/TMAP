/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_DEBUG_EXACT_H_
#define TMAP_DEBUG_EXACT_H_

/*! 
  Debugging Functions
  */ 

/*! 
  structure to store the command line options for 'tmap exact'
  */
typedef struct {
    char **argv;  /*!< the command line argv structure */
    int argc;  /*!< the number of command line arguments passed */
    char *fn_fasta;  /*!< the fasta reference file name (-f) */
    char *fn_reads;  /*!< the fastq reads file name (-r) */
    int32_t n_only;  /*!< print the number of matches only (-c) */
} tmap_debug_exact_opt_t;
/*! 
  main-like function for 'tmap exact'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int tmap_debug_exact(int argc, char *argv[]);
#endif
