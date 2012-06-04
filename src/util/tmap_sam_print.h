/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_SAM_PRINT_H
#define TMAP_SAM_PRINT_H

#include <config.h>
#include "../samtools/bam.h"
#include "../index/tmap_refseq.h"
#include "../io/tmap_file.h"
#include "../io/tmap_seq_io.h"

/*! 
*/
                    
#define TMAP_SAM_PRINT_VERSION "1.4"

/*! 
  prints out a SAM header
  @param  fp            the output file pointer
  @param  refseq        pointer to the reference sequence (forward)
  @param  seqio         the input reading data structure, NULL otherwise
  @param  sam_rg        the SAM RG line, NULL otherwise
  @param  sam_flowspace_tags  1 if SFF specific SAM tags are to be outputted, 0 otherwise
  @param  ignore_rg_sam_tags 1 if we are to ignore RG tags from the input file, 0 otherwise
  @param  argc          the number of input command line arguments
  @param  argv          the input command line arguments
  @details              the following header tags will be outptted: \@SQ:SN:LN and \@PG:ID:VN:CL.  
  The sam_flowspace_tags tags precedence over ignore_rg_sam_tags
  */
void
tmap_sam_print_header(tmap_file_t *fp, tmap_refseq_t *refseq, 
                      tmap_seq_io_t *seqio, char *sam_rg, 
                      int32_t sam_flowspace_tags, int32_t ignore_rg_sam_tags, 
                      int argc, char *argv[]);

/*! 
  prints out a SAM record signifying the sequence is unmapped 
  @param  fp            the file pointer to which to print
  @param  seq           the sequence that is unmapped
  @param  sam_flowspace_tags  1 if SFF specific SAM tags are to be outputted, 0 otherwise
  @param  bidirectional  1 if a bidirectional SAM tag is to be added, 0 otherwise
  @param  refseq        pointer to the reference sequence (forward)
  @param  end_num       0 if there is no mate (all mate params ignored), 1 if the mate is the first fragment, 2 if the mate is the last fragment
  @param  m_unmapped    1 if the mate is unmapped, 0 otherwise (m_strand/m_seqid/m_pos are ignored)
  @param  m_prop        1 if properly paired, 0 otherwise
  @param  m_strand      the mates strand
  @param  m_seqid       the mates seqid (zero-based), 0 otherwise
  @param  m_pos         the mates position (zero-based), 0 otherwise
  @param  format      optional tag format (printf-style)
  @param  ...         arguments for the format
  */
inline void
tmap_sam_print_unmapped(tmap_file_t *fp, tmap_seq_t *seq, int32_t sam_flowspace_tags, int32_t bidirectional, tmap_refseq_t *refseq,
                        uint32_t end_num, uint32_t m_unmapped, uint32_t m_prop, 
                        uint32_t m_strand, uint32_t m_seqid, uint32_t m_pos,
                        const char *format, ...);


/*! 
  prints out a mapped SAM record 
  @param  fp          the file pointer to which to print
  @param  seq         the sequence that is mapped
  @param  sam_flowspace_tags  1 if SFF specific SAM tags are to be outputted, 0 otherwise
  @param  bidirectional  1 if a bidirectional SAM tag is to be added, 0 otherwise
  @param  seq_eq      1 if the SEQ field is to use '=' symbols, 0 otherwise
  @param  refseq      pointer to the reference sequence (forward)
  @param  strand      the strand of the mapping
  @param  seqid       the sequence index (0-based)
  @param  pos         the position (0-based)
  @param  secondary   1 if the alignment is a secondary alignment, 0 otherwise
  @param  end_num     0 if there is no mate (all mate params ignored), 1 if the mate is the first fragment, 2 if the mate is the last fragment
  @param  m_unmapped  1 if the mate is unmapped, 0 otherwise (m_strand/m_seqid/m_pos/m_tlen are ignored)
  @param  m_prop      1 if properly paired, 0 otherwise
  @param  m_num_std   the number of standard devaitions from the mean insert size if paired
  @param  m_strand    the mates strand
  @param  m_seqid     the mates seqid (zero-based), 0 otherwise
  @param  m_pos       the mates position (zero-based), 0 otherwise
  @param  m_tlen      the mate template length (zero-based), 0 otherwise
  @param  mapq        the mapping quality
  @param  cigar       the cigar array
  @param  n_cigar     the number of cigar operations
  @param  score       the alignment score
  @param  ascore      the original base alignment score (SFF only)
  @param  pscore      the pairing alignment score (paired reads only)
  @param  nh          the number of reported alignments (NH tag)
  @param  algo_id     the algorithm id
  @param  algo_stage  the algorithm stage (1 or 2) 
  @param  format      optional tag format (printf-style)
  @param  ...         arguments for the format
  @details            the format should not include the MD tag, which will be outputted automatically
  */
inline void
tmap_sam_print_mapped(tmap_file_t *fp, tmap_seq_t *seq, int32_t sam_flowspace_tags, int32_t bidirectional, int32_t seq_eq, tmap_refseq_t *refseq,
                      uint8_t strand, uint32_t seqid, uint32_t pos, int32_t secondary,
                      uint32_t end_num, uint32_t m_unmapped, uint32_t m_prop, double m_num_std, uint32_t m_strand,
                      uint32_t m_seqid, uint32_t m_pos, uint32_t m_tlen,
                      uint8_t mapq, uint32_t *cigar, int32_t n_cigar,
                      int32_t score, int32_t ascore, int32_t pscore, int32_t nh, int32_t algo_id, int32_t algo_stage,
                      const char *format, ...);

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
