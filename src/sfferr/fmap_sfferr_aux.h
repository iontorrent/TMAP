#ifndef FMAP_SFFERR_AUX_H_
#define FMAP_SFFERR_AUX_H_

#include <config.h>
#include "../index/fmap_refseq.h"
#include "../seq/fmap_sff.h"
#ifdef HAVE_SAMTOOLS
#include <sam.h>
#endif

#ifdef HAVE_SAMTOOLS

/*!
 Stores the empirical probibility distribution: Pr(flow signal | polymer call, flow index)
 */
typedef struct {
    uint16_t ***signal_counts;  /*!< [i][j][k] => [flow_length][max_polymer][max_signal] */
    int32_t flow_length;  /*!< the maximum number of flows */
    int32_t max_polymer;  /*!< the maximum length called polymer */
    int32_t max_signal;  /*!< the maximum valued signal (2^16)-1 */
} fmap_sfferr_aux_pr_flow_given_call_t;

/*!
  initialization structure
  @param  header  pointer the SFF header
  @return         pointer to the initialized structure
  @details        the header gives the flow length
  */
fmap_sfferr_aux_pr_flow_given_call_t *
fmap_sfferr_aux_pr_flow_given_call_init(fmap_sff_header_t *header);

/*!
  destroy
  @param  ptr  a pointer to the structure
  */
void
fmap_sfferr_aux_pr_flow_given_call_destroy(fmap_sfferr_aux_pr_flow_given_call_t *ptr);

/*!
  reallocates the given structure 
  @param  ptr  a pointer to the structure
  @param  flow_length   the maximum flow length
  @param  max_polymer  the maximum polymer length
  @details     this will only expand the structure
  */
void
fmap_sfferr_aux_pr_flow_given_call_realloc(fmap_sfferr_aux_pr_flow_given_call_t *ptr, int32_t flow_length, int32_t max_polymer);

/*!
  reads the structure from disk
  @param  fp  a file pointer from which to read
  @return     a pointer to the structure
  */
fmap_sfferr_aux_pr_flow_given_call_t *
fmap_sfferr_aux_pr_flow_given_call_read(fmap_file_t *fp);

/*!
  write the structure to disk
  @param  fp   a file pointer to which to write
  @param  ptr  a pointer to the structure to write
  */
void
fmap_sfferr_aux_pr_flow_given_call_write(fmap_sfferr_aux_pr_flow_given_call_t *ptr, fmap_file_t *fp);

/*!
  updates the given empirical probability distributions
  @param  refseq  pointer to the reference sequence structure
  @param  sff     pointer to the sff structure
  @param  bam     pointer to the bam structure
  @param  ptr1    pointer to the first probability distribution
  */
void
fmap_sfferr_aux(fmap_refseq_t *refseq, fmap_sff_t *sff, bam1_t *bam,
                fmap_sfferr_aux_pr_flow_given_call_t *ptr1);

#endif /* HAVE_SAMTOOLS */

#endif
