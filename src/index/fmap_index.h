#ifndef FMAP_INDEX_H_
#define FMAP_INDEX_H_

#define FMAP_INDEX_LARGE_GENOME 0x40000000 
#define FMAP_INDEX_TOO_BIG_GENOME 0x100000000 

/*! 
  details  Constructs the packed reference sequence, BWT string, and Suffix Array.
  */

/*! 
            structure to store the command line options for 'fmap index'
  @field  fn_fasta      the fasta file name (-f)
  @field  occ_interval  the occurrence array interval (-o)
  @field  hash_width    the occurrence hash width (-w)
  @field  sa_interval   the suffix array interval (-i)
  @field  is_large      0 to use the short BWT construction algorith, 1 otherwise (large BWT construction algorithm)
  */
typedef struct {
    char *fn_fasta;
    int32_t occ_interval;
    int32_t hash_width;
    int32_t sa_interval;
    int32_t is_large;
} fmap_index_opt_t;

/*! 
     main-like function for 'fmap index'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int fmap_index(int argc, char *argv[]);

#endif
