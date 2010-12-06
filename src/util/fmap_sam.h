#ifndef FMAP_SAM_H_
#define FMAP_SAM_H_

#include <config.h>
#ifdef HAVE_SAMTOOLS
#include <bam.h>
#endif
#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../io/fmap_file.h"

/*! 
*/

/*! 
  prints out a SAM header
  @param  fp      the output file pointer
  @param  refseq  pointer to the reference sequence (forward)
  @param  argc    the number of input command line arguments
  @param  argv    the input command line arguments
  @details        the following header tags will be ouptted: \@SQ:SN:LN and \@PG:ID:VN:CL.
  */
void
fmap_sam_print_header(fmap_file_t *fp, fmap_refseq_t *refseq, int argc, char *argv[]);

/*! 
  prints out a SAM record signifying the sequence is unmapped 
  @param  fp   the file pointer to which to print
  @param  seq  the sequence that is unmapped
  */
inline void
fmap_sam_print_unmapped(fmap_file_t *fp, fmap_seq_t *seq);

/*! 
  prints out a mapped SAM record 
  @param  fp       the file pointer to which to print
  @param  seq      the sequence that is mapped
  @param  refseq   pointer to the reference sequence (forward)
  @param  strand   the strand of the mapping
  @param  seqid    the sequence index (0-based)
  @param  pos      the position (0-based)
  @param  mapq     the mapping quality
  @param  cigar    the cigar array
  @param  n_cigar  the number of cigar operations
  @param  format   optional tag format (printf-style)
  @param  ...      arguments for the format
  @details         the format should not include the MD tag, which will be outputted automatically
  */
inline void
fmap_sam_print_mapped(fmap_file_t *fp, fmap_seq_t *seq, fmap_refseq_t *refseq,
                      uint8_t strand, uint32_t seqid, uint32_t pos,
                      uint8_t mapq, uint32_t *cigar, int32_t n_cigar,
                      const char *format, ...);

#ifdef HAVE_SAMTOOLS
/*!
  recreates an MD given the new reference/read alignment
  @param  b     the SAM/BAM structure
  @param  ref   the reference
  @param  len   the length of the alignment
  */
void 
fmap_sam_md1(bam1_t *b, char *ref, int32_t len);

/*!
  updates the cigar and MD given the new reference/read alignment
  @param  b     the SAM/BAM structure
  @param  ref   the reference
  @param  read  the read
  @param  len   the length of the alignment
  */
void
fmap_sam_update_cigar_and_md(bam1_t *b, char *ref, char *read, int32_t len);
#endif

#endif
