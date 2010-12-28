/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SAM_H
#define TMAP_SAM_H

#include <config.h>
#ifdef HAVE_SAMTOOLS
#include <bam.h>
#endif
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../io/tmap_file.h"
#include "../io/tmap_seq_io.h"

/*! 
*/
                    
#define TMAP_SAM_VERSION "1.3"

/*! 
  prints out a SAM header
  @param  fp            the output file pointer
  @param  refseq        pointer to the reference sequence (forward)
  @param  seqio         the input reading data structure, NULL otherwise
  @param  sam_rg        the SAM RG line, NULL otherwise
  @param  sam_sff_tags  1 if SFF specific SAM tags are to be outputted, 0 otherwise
  @param  argc          the number of input command line arguments
  @param  argv          the input command line arguments
  @details              the following header tags will be ouptted: \@SQ:SN:LN and \@PG:ID:VN:CL.
  */
void
tmap_sam_print_header(tmap_file_t *fp, tmap_refseq_t *refseq, 
                      tmap_seq_io_t *seqio, char *sam_rg, int32_t sam_sff_tags, 
                      int argc, char *argv[]);

/*! 
  prints out a SAM record signifying the sequence is unmapped 
  @param  fp            the file pointer to which to print
  @param  seq           the sequence that is unmapped
  @param  sam_sff_tags  1 if SFF specific SAM tags are to be outputted, 0 otherwise
  */
inline void
tmap_sam_print_unmapped(tmap_file_t *fp, tmap_seq_t *seq, int32_t sam_sff_tags);

/*! 
  prints out a mapped SAM record 
  @param  fp          the file pointer to which to print
  @param  seq         the sequence that is mapped
  @param  sam_sff_tags  1 if SFF specific SAM tags are to be outputted, 0 otherwise
  @param  refseq      pointer to the reference sequence (forward)
  @param  strand      the strand of the mapping
  @param  seqid       the sequence index (0-based)
  @param  pos         the position (0-based)
  @param  mapq        the mapping quality
  @param  cigar       the cigar array
  @param  n_cigar     the number of cigar operations
  @param  score       the alignment score
  @param  ascore      the original base alignment score (SFF only)
  @param  algo_id     the algorithm id
  @param  algo_stage  the algorithm stage (1 or 2) 
  @param  format      optional tag format (printf-style)
  @param  ...         arguments for the format
  @details            the format should not include the MD tag, which will be outputted automatically
  */
inline void
tmap_sam_print_mapped(tmap_file_t *fp, tmap_seq_t *seq, int32_t sam_sff_tags, tmap_refseq_t *refseq,
                      uint8_t strand, uint32_t seqid, uint32_t pos,
                      uint8_t mapq, uint32_t *cigar, int32_t n_cigar,
                      int32_t score, int32_t ascore, int32_t algo_id, int32_t algo_stage,
                      const char *format, ...);

#ifdef HAVE_SAMTOOLS
/*!
  recreates an MD given the new reference/read alignment
  @param  b     the SAM/BAM structure
  @param  ref   the reference
  @param  len   the length of the alignment
  */
void 
tmap_sam_md1(bam1_t *b, char *ref, int32_t len);

/*!
  updates the cigar and MD given the new reference/read alignment
  @param  b     the SAM/BAM structure
  @param  ref   the reference
  @param  read  the read
  @param  len   the length of the alignment
  */
void
tmap_sam_update_cigar_and_md(bam1_t *b, char *ref, char *read, int32_t len);
#endif

#endif
