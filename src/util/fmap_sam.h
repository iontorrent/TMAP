#ifndef FMAP_SAM_H_
#define FMAP_SAM_H_

#include "../seq/fmap_seq.h"
#include "../index/fmap_refseq.h"
#include "../io/fmap_file.h"

/*! @function
  @abstract  prints out a SAM header
  @param  fp      the output file pointer
  @param  refseq  pointer to the reference sequence (forward)
  @param  argc    the number of input command line arguments
  @param  argc    the input command line arguments
  @discussion     the following header tags will be ouptted: @SQ:SN:LN and @PG:ID:VN:CL.
  */
void
fmap_sam_print_header(fmap_file_t *fp, fmap_refseq_t *refseq, int argc, char *argv[]);

/*! @function
  @abstract  prints out a SAM record signifying the sequence is unmapped 
  @param  fp   the file pointer to which to print
  @param  seq  the sequence that is unmapped
  */
inline void
fmap_sam_print_unmapped(fmap_file_t *fp, fmap_seq_t *seq);

#endif
