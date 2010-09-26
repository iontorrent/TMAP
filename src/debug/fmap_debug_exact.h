#ifndef FMAP_DEBUG_EXACT_H_
#define FMAP_DEBUG_EXACT_H_

/*! 
  Debugging Functions
 */ 

/*! 
  structure to store the command line options for 'fmap exact'
*/
typedef struct {
  char *fn_fasta;  /*!< the fasta reference file name (-f) */
  char *fn_reads;  /*!< the fastq reads file name (-r) */
} fmap_debug_exact_opt_t;
/*! 
  main-like function for 'fmap exact'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
*/
int fmap_debug_exact(int argc, char *argv[]);
#endif
