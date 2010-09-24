#ifndef FMAP_MAP2_H_
#define FMAP_MAP2_H_

#include <stdlib.h>
#include <sys/ipc.h>

/*! @typedef
  @abstract
  @abstract                structure to store the command line options for 'fmap map2'
  @field  argv              the command line argv structure
  @field  argc              the number of command line arguments passed
  @field  fn_fasta          the fasta reference file name (-f)
  @field  fn_reads          the reads file name (-r)
  @field  reads_format      the reads file format (-F)
  @field  score_match         the match score (-A)
  @field  pen_mm            the mismatch penalty (-M)
  @field  pen_gapo          the indel open penalty (-O)
  @field  pen_gape          the indel extension penalty (-E)
  @field  yita              the error recurrence coefficient (-y) 
  @field  mask_level        the mask level (-m)
  @field  length_coef       the coefficient of length-threshold adjustment (-c)
  @field  band_width        the band width (-w) 
  @field  score_thr         the score threshold (match-score-scaled) (-T)
  @field  max_seed_intv     the maximum seed interval (-S)
  @field  z_best            the number of top scoring hits to keep (-Z)
  @field  seeds_rev         the maximum number of seeds for which reverse alignment is triggered (-N)
  @field  reads_queue_size  the reads queue size (-q)
  @field  num_threads       the number of threads (-n)
  @field  input_compr       the input compression type (-j and -z)
  @field  output_compr      the output compression type (-J and -Z)
  @field  shm_key           the shared memory key (-s)
  */
typedef struct {
    char **argv;
    int argc;
    char *fn_fasta;
    char *fn_reads;
    int32_t reads_format;
    int32_t score_match, pen_mm, pen_gapo, pen_gape;
    double yita, mask_level, length_coef;
    int32_t band_width, score_thr, max_seed_intv, z_best, seeds_rev;
    int32_t reads_queue_size, num_threads, input_compr, output_compr;
    key_t shm_key;
} fmap_map2_opt_t;

/*! @typedef
  @abstract                 data to be passed to a thread
  @field  seq_buffer         the buffer of sequences
  @field  seq_buffer_length  the buffer length
  @field  sams               the sam alignments for each sequence
  @field  refseq             pointer to the reference sequence (forward)
  @field  bwt                pointer to the BWT indices (forward/reverse)
  @field  sa                 pointer to the SA (forward/reverse)
  @field  tid                the zero-based thread id
  @field  opt                the options to this program
  */
typedef struct {
    fmap_seq_t **seq_buffer;
    int32_t seq_buffer_length;
    struct __fmap_map2_sam_t **sams;
    fmap_refseq_t *refseq;
    fmap_bwt_t *bwt[2];
    fmap_sa_t *sa[2];
    int32_t tid;
    fmap_map2_opt_t *opt;
} fmap_map2_thread_data_t;

/*! @function
  @abstract     main-like function for 'fmap map2'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
fmap_map2_main(int argc, char *argv[]);

#endif
