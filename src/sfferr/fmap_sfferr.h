#ifndef FMAP_SFFERR_H_
#define FMAP_SFFERR_H_

#include <config.h>

#ifdef HAVE_SAMTOOLS

/*! 
  @brief                structure to store the command line options for 'fmap sfferr'
  */
typedef struct {
    char **argv;  /*!< the command line argv structure */
    int argc;  /*!< the number of command line arguments passed */
    char *fn_fasta;  /*!< the fasta reference file name (-f) */
    char *fn_sff;  /*!< the sff file name (-r) */
    char *fn_sam;  /*!< the SAM file name (-S) */
    int32_t reads_format;  /*!< the reads file format (-F) */
    double rand_sample_num;  /*!< the fraction of reads to sample (-R) */
    int32_t input_compr;  /*!< the input compression type (-j and -z) */
    key_t shm_key;  /*!< the shared memory key (-s) */
} fmap_sfferr_opt_t;

/*! 
  @brief     main-like function for 'fmap sfferr'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int fmap_sfferr_main(int argc, char *argv[]);
#endif /* HAVE_SAMTOOLS */

#endif
