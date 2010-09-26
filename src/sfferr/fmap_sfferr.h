#ifndef FMAP_SFFERR_H_
#define FMAP_SFFERR_H_

/*! 
  @brief                structure to store the command line options for 'fmap sfferr'
*/
typedef struct {
  char **argv;  /*! the command line argv structure */
  int argc;  /*! the number of command line arguments passed */
  char *fn_reads;  /*! the reads file name (-r) */
  int32_t reads_format;  /*! the reads file format (-F)  */
  double rand_sample_num;  /*! the fraction of reads to sample (-R) */
  int32_t input_compr;  /*! the input compression type (-j and -z) */
} fmap_sff_err_opt_t;
/*! 
  @brief     main-like function for 'fmap sfferr'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
*/
int fmap_map1(int argc, char *argv[]);
#endif
