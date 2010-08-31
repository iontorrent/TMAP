#ifndef FMAP_DEBUG_HASH_H_
#define FMAP_DEBUG_HASH_H_

/*! @typedef
  @abstract        structure to store the command line options for 'fmap hash'
  @field  fn_fasta  the fasta reference file name (-f)
*/
typedef struct {
    char *fn_fasta;
    int32_t hash_width;
} fmap_debug_hash_opt_t;

/*! @function
  @abstract     main-like function for 'fmap hash'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
*/
int fmap_debug_hash(int argc, char *argv[]);

#endif
