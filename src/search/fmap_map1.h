#ifndef FMAP_MAP1_H_
#define FMAP_MAP1_H_

/*! @typedef
  @abstract                structure to store the command line options for 'fmap exact'
  @field  fn_fasta          the fasta reference file name (-f)
  @field  fn_reads          the reads file name (-r)
  @field  reads_format      the reads file format (-F) 
  @field  seed_length       the k-mer length to seed CALs (-s)
  @field  max_mm            maximum number of mismatches (-m)
  @field  max_mm_frac       maximum (read length) fraction of mismatches (-m)
  @field  max_gapo          maximum number of indel opens (-o)
  @field  max_gapo_frac     maximum (read length) fraction of indel opens (-o)
  @field  max_gape          maximum number of indel extensions (-e)
  @field  max_gape_frac     maximum fraction of indel extensions (-e)
  @field  pen_mm            the mismatch penalty (-M)
  @field  pen_gapo          the indel open penalty (-O)
  @field  pen_gape          the indel extension penalty (-E)
  @field  max_cals_del      the maximum number of CALs to extend a deletion (-d)
  @field  indel_ends_bound  indels are not allowed within INT number of bps from the end of the read (-i)
  @field  max_best_cals     stop searching when INT optimal CALs have been found (-b)
  @field  reads_queue_size  the reads queue size (-q)
  @field  num_threads       the number of threads (-n)
*/
typedef struct {
    char *fn_fasta;
    char *fn_reads;
    int32_t reads_format;
    int32_t seed_length;
    int32_t max_mm;
    double max_mm_frac;
    int32_t max_gapo;
    double max_gapo_frac;
    int32_t max_gape;
    double max_gape_frac;
    int32_t pen_mm;
    int32_t pen_gapo;
    int32_t pen_gape;
    int32_t max_cals_del;
    int32_t indel_ends_bound;
    int32_t max_best_cals;
    int32_t reads_queue_size;
    int32_t num_threads;
} fmap_map1_opt_t;

/*! @typedef
  @abstract                 data to be passed to a thread
  @field  seq_buffer         the buffer of sequences
  @field  seq_buffer_length  the buffer length
  @field  refseq             pointer to the references sequence
  @field  bwt                pointer to the BWT index
  @field  sa                 pointer to the SA structure
  @field  tid                the zero-based thread id
  @field  opt                the options to this program
 */
typedef struct {
    fmap_seq_t **seq_buffer;
    int32_t seq_buffer_length;
    fmap_refseq_t *refseq;
    fmap_bwt_t *bwt;
    fmap_sa_t *sa;
    int32_t tid;
    fmap_map1_opt_t *opt;
} fmap_map1_thread_data_t;

/*! @function
  @abstract     main-like function for 'fmap exact'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
*/
int fmap_map1(int argc, char *argv[]);

#endif
