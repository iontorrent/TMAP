/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sam_print.h"
#include "../samtools/sam_header.h"
#include "tmap_seq_io.h"
#include "tmap_sff_io.h"
#include "tmap_seq_io.h"
#include "tmap_seqs_io.h"

inline tmap_seqs_io_t*
tmap_seqs_io_init(char **fns, int32_t fn_num, int8_t seq_type, int32_t compression)
{
  tmap_seqs_io_t *io= NULL;
  int32_t i;

  io = tmap_calloc(1, sizeof(tmap_seqs_io_t), "io");
  io->type = seq_type;
      
  if(1 < io->n && (TMAP_SEQ_TYPE_SAM == io->type || TMAP_SEQ_TYPE_BAM == io->type)) {
      tmap_error("Multi-SAM/BAM not supported", Exit, OutOfRange);
  }

  if(NULL == fns) { // stdin
      io->n = 1;
      io->seqios = tmap_calloc(1, sizeof(tmap_seq_io_t*), "io->seqios");
      io->seqios[0] = tmap_seq_io_init("-", seq_type, 0, compression); // NB: always reading
  }
  else { // from file(s)
      io->n = fn_num;
      io->seqios = tmap_calloc(fn_num, sizeof(tmap_seq_io_t*), "io->seqios");
      for(i=0;i<io->n;i++) {
          io->seqios[i] = tmap_seq_io_init(fns[i], seq_type, 0, compression); // NB: always reading
      }
  }

  return io;
}

inline void
tmap_seqs_io_destroy(tmap_seqs_io_t *io)
{
  int32_t i;
  for(i=0;i<io->n;i++) {
      tmap_seq_io_destroy(io->seqios[i]);
  }
  free(io->seqios);
  free(io);
}

inline int
tmap_seqs_io_read(tmap_seqs_io_t *io, tmap_seqs_t *seqs)
{
  int32_t i;

  /*
   * Case 1 - SAM/BAM
   *    - NB: there must only be one input file
   *    - Read a record, if paired, then read the next
   * Case 2 - SFF/FQ
   *    - NB: there can be zero or more input files
   *    - Read one from each file, store in one record
   */

  if(io->type != seqs->type) {
      tmap_error("type mismatch", Exit, OutOfRange);
  }

  // reset seqs
  seqs->n = 0;
  if(TMAP_SEQ_TYPE_SAM == io->type || TMAP_SEQ_TYPE_BAM == io->type) {
      // NB: to supported paired reads, we check the paired flag
      for(i=0;i<2;i++) {
          tmap_seq_t *seq = tmap_seqs_get(seqs, i);
          if(tmap_seq_io_read(io->seqios[0], seq) < 0) return EOF; // TODO: better error checking
          tmap_seqs_add(seqs, seq); 
          // break if not paired
          if(0 == (seq->data.sam->b->core.flag & BAM_FPAIRED)) break;
      }
  }
  else {
      // read in one per file
      for(i=0;i<io->n;i++) {
          tmap_seq_t *seq = tmap_seqs_get(seqs, i);
          if(tmap_seq_io_read(io->seqios[i], seq) < 0) return EOF; // TODO: better error checking
          tmap_seqs_add(seqs, seq); 
      }
  }

  return 0;
}

int
tmap_seqs_io_read_buffer(tmap_seqs_io_t *io, tmap_seqs_t **seqs_buffer, int32_t buffer_length)
{
  int32_t n = 0;

  if(buffer_length <= 0) return 0;

  while(n < buffer_length) {
      if(NULL == seqs_buffer[n]) {
          seqs_buffer[n] = tmap_seqs_init(io->type);
      }
      else if(io->type != seqs_buffer[n]->type) { // check if we need to change the type
          tmap_seqs_destroy(seqs_buffer[n]);
          seqs_buffer[n] = tmap_seqs_init(io->type);
      }
      if(tmap_seqs_io_read(io, seqs_buffer[n]) < 0) {
          break;
      }
      n++;
  }

  return n;
}

static void
tmap_seqs_io_init2_fs_and_add(tmap_seqs_io_t *io_in,
                         sam_header_records_t *records,
                         sam_header_record_t *record,
                         int32_t sam_flowspace_tags)
{
  char tag[2];
  // add @RG.KS and @RG.FO
  if(io_in->type == TMAP_SEQ_TYPE_SFF && sam_flowspace_tags) {
      if(io_in->n <= records->n) tmap_error("Too many read groups specified", Exit, OutOfRange);
      // @RG.KS
      tag[0]='K';tag[1]='S';
      if(0 == sam_header_record_add(record, tag, tmap_sff_io_get_rg_ks(io_in->seqios[records->n]->io.sffio))) {
          tmap_error("Could not add the KS tag", Exit, OutOfRange);
      }
      // @RG.FO
      tag[0]='F';tag[1]='O';
      if(0 == sam_header_record_add(record, tag, tmap_sff_io_get_rg_fo(io_in->seqios[records->n]->io.sffio))) {
          tmap_error("Could not add the FO tag", Exit, OutOfRange);
      }
  }
  // check for the @RG.ID and @RG.SM tags
  if(NULL == sam_header_record_get(record, "ID")) tmap_bug(); // should not happen
  if(NULL == sam_header_record_get(record, "SM")) tmap_error("SM tag missing from the read group", Exit, OutOfRange);
  // add the read group
  if(0 == sam_header_records_add(records, record)) tmap_bug(); 
}

bam_header_t *
tmap_seqs_io_to_bam_header(tmap_refseq_t *refseq,
                           tmap_seqs_io_t *io_in,
                           char **rg_sam, int32_t rg_sam_num,
                           int32_t sam_flowspace_tags,
                           int32_t argc, char *argv[])
{
  bam_header_t *bam_header = NULL;
  sam_header_t *header = NULL; // the output header
  sam_header_records_t *records = NULL;
  sam_header_record_t *record = NULL;
  char tag[2];
  char *command_line= NULL;
  int32_t i, j;

  // @HD
  if(io_in->type == TMAP_SEQ_TYPE_SAM || io_in->type == TMAP_SEQ_TYPE_BAM) {
      // should be only one input file
      if(1 != io_in->n) {
          tmap_bug();
      }
      // get the current header
      header = io_in->seqios[0]->io.samio->fp->header->header; // wow, that's a lot of pointers
      if(NULL == header) {
          tmap_bug();
      }
  }
  else {
      // empty header
      header = sam_header_init();
      // @HD - header line
      records = sam_header_get_records(header, "HD"); // get the header line
      if(NULL == records) tmap_bug();
      record = sam_header_record_init("HD"); // new header line
      if(0 == sam_header_record_add(record, "VN", "1.4")) tmap_bug(); // version number
      if(0 == sam_header_records_add(records, record)) tmap_bug(); // add the header line
      // nullify
      record = NULL;
      records = NULL;
  }

  // @SQ
  if(NULL != refseq) {
      records = sam_header_get_records(header, "SQ"); // get the sequence dictionary
      if(NULL == records) tmap_bug();
      for(i=0;i<refseq->num_annos;i++) { // for each reference sequence
          char num[32];
          record = sam_header_record_init("SQ"); // new reference sequence record
          if(0 == sam_header_record_add(record, "SN", refseq->annos[i].name->s)) tmap_bug(); // reference sequence name
          if(sprintf(num, "%u", (uint32_t)refseq->annos[i].len) < 0) tmap_bug(); // integer to string
          if(0 == sam_header_record_add(record, "SN", num)) tmap_bug(); // reference sequence name
          if(0 == sam_header_records_add(records, record)) tmap_bug(); // add the reference sequence record
      }
      record = NULL;
      records = NULL;
  }

  // @RG - read group
  if(0 < rg_sam_num) { // @RG specified on the command line
      // Check for SAM/BAM
      if(io_in->type == TMAP_SEQ_TYPE_SAM || io_in->type == TMAP_SEQ_TYPE_BAM) {
          tmap_error("Cannot specify the read groups on the command line when using SAM/BAM as input."
                     "  Please embed in the SAM/BAM header instead.", Exit, OutOfRange);
      }
      record = NULL;
      records = sam_header_get_records(header, "RG"); // get the read group
      if(NULL == records) tmap_bug();
      // go through the command line arguments
      for(i=0;i<rg_sam_num;i++) {
          if(strlen(rg_sam[i]) < 4) tmap_error("Read group too small", Exit, OutOfRange);
          if(':' != rg_sam[i][2]) tmap_error("Read group improperly formatted (no colon)", Exit, OutOfRange);

          // check for id
          if('I' == rg_sam[i][0] && 'D' == rg_sam[i][1]) { // new read group
              if(NULL != record) { // add the record
                  tmap_seqs_io_init2_fs_and_add(io_in, records, record, sam_flowspace_tags); // add @RG.KS and @RG.FO
              }
              record = sam_header_record_init("RG"); // new read group
          }
          // add the tag/value to the record
          tag[0]=rg_sam[i][0]; tag[1]=rg_sam[i][1]; // setup the tag
          if(0 == sam_header_record_add(record, tag, rg_sam[i]+3)) tmap_bug(); // add the tag/value
      }
      if(NULL != record) { // add the record
          tmap_seqs_io_init2_fs_and_add(io_in, records, record, sam_flowspace_tags); // add @RG.KS and @RG.FO
      }
      // check that the # of read groups added was the same as the # of input files...
      if(records->n != io_in->n) tmap_error("The number of read groups did not match the number of input files", Exit, OutOfRange);
  }
  else if(io_in->type != TMAP_SEQ_TYPE_SAM && io_in->type != TMAP_SEQ_TYPE_BAM) { // dummy...
      records = sam_header_get_records(header, "RG"); // get the read group
      if(NULL == records) tmap_bug();
      for(i=0;i<io_in->n;i++) { // for each input file
          char buf[32];
          record = sam_header_record_init("RG"); // new read group
          if(1 == io_in->n) strcpy(buf, "NOID");
          else if(sprintf(buf, "NOID.%d", i+1) < 0) tmap_bug();
          if(0 == sam_header_record_add(record, "ID", buf)) tmap_bug(); // dummy ID
          if(0 == sam_header_record_add(record, "SM", "NOSM")) tmap_bug(); // dummy SM, for Picard validation
          tmap_seqs_io_init2_fs_and_add(io_in, records, record, sam_flowspace_tags); // add @RG.KS and @RG.FO
      }
  }

  // @PG - program group
  // TODO: check for previous program group ID and set @PG.PP
  records = sam_header_get_records(header, "PG"); // get the program group records
  if(NULL == records) tmap_bug();
  record = sam_header_record_init("PG"); // new program group
  if(0 == sam_header_record_add(record, "ID", PACKAGE_NAME)) tmap_bug(); // @PG.ID
  if(0 == sam_header_record_add(record, "VN", PACKAGE_VERSION)) tmap_bug(); // @PG.VN
  // @PG.CL
  command_line = NULL;
  j = 1; // for the EOL
  for(i=0;i<argc;i++) {
      if(0 < i) j++;
      j += strlen(argv[i]);
      command_line = tmap_realloc(command_line, sizeof(char) * j, "command_line");
      if(0 < i) strcat(command_line, " ");
      strcat(command_line, argv[i]);
  }
  if(0 == sam_header_record_add(record, "CL", command_line)) tmap_bug(); // @PG.CL
  free(command_line);


  // Create a BAM Header from the SAM Header
  bam_header = bam_header_init(); // empty
  bam_header->header = header; // soft-copy the header
  bam_header = sam_header_to_bam_header(bam_header); // convert

  return bam_header;
}
