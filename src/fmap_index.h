#ifndef FMAP_INDEX_H_
#define FMAP_INDEX_H_

/*! @typedef
  @abstract        structure to store the command line options for 'fmap index'
  @field  fn_fasta  the fasta file name (-f)
*/
typedef struct {
    char *fn_fasta;
} fmap_index_opt_t;

/*! @function
  @abstract     main-like function for 'fmap index'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
*/
int fmap_index(int argc, char *argv[]);

#endif
