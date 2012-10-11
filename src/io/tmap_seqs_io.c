/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sam_convert.h"
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
tmap_seqs_io_read(tmap_seqs_io_t *io, tmap_seqs_t *seqs, sam_header_t *header)
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
          tmap_seq_update(seq, i, header);
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
          tmap_seq_update(seq, i, header);
      }
  }

  return 0;
}

int
tmap_seqs_io_read_buffer(tmap_seqs_io_t *io, tmap_seqs_t **seqs_buffer, int32_t buffer_length, sam_header_t *header)
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
      if(tmap_seqs_io_read(io, seqs_buffer[n], header) < 0) {
          break;
      }
      n++;
  }

  return n;
}

static void
tmap_seqs_io_init2_fs_and_add(tmap_seqs_io_t *io_in,
                              sam_header_t *header,
                              sam_header_record_t *record)
{
  char tag[2];
  // add @RG.KS and @RG.FO
  if(io_in->type == TMAP_SEQ_TYPE_SFF) {
      sam_header_records_t *records = sam_header_get_records(header, record->tag); // get the header line
      if(io_in->n <= records->n) tmap_error("Too many read groups specified", Exit, OutOfRange);
      // @RG.KS
      tag[0]='K';tag[1]='S';
      if(0 == sam_header_record_add(record, tag, tmap_sff_io_get_rg_ks(io_in->seqios[records->n]->io.sffio))) {
          tmap_error("Could not add the KS tag; most likely it is already present", Exit, OutOfRange);
      }
      // @RG.FO
      tag[0]='F';tag[1]='O';
      if(0 == sam_header_record_add(record, tag, tmap_sff_io_get_rg_fo(io_in->seqios[records->n]->io.sffio))) {
          tmap_error("Could not add the FO tag; most likely it is already present", Exit, OutOfRange);
      }
  }
  // check for the @RG.ID and @RG.SM tags
  if(NULL == sam_header_record_get(record, "ID")) tmap_bug(); // should not happen
  if(NULL == sam_header_record_get(record, "SM")) {
      if(0 == sam_header_record_add(record, "SM", "NOSM")) tmap_bug(); // dummy SM, for Picard validation
  }
  if(NULL == sam_header_record_get(record, "PG")) {
      if(0 == sam_header_record_add(record, "PG", PACKAGE_NAME)) tmap_bug(); // dummy PG
  }
  // add the read group
  if(0 == sam_header_add_record(header, record)) tmap_bug(); 
}

bam_header_t *
tmap_seqs_io_to_bam_header(tmap_refseq_t *refseq,
                           tmap_seqs_io_t *io_in,
                           char **rg_sam, int32_t rg_sam_num,
                           int32_t argc, char *argv[])
{
  bam_header_t *bam_header = NULL;
  sam_header_t *header = NULL; // the output header
  sam_header_record_t *record = NULL;
  sam_header_record_t **record_list = NULL;
  char tag[2];
  char *command_line= NULL;
  char *id = NULL;
  char *id_pp = NULL;
  int32_t i, j;

  // @HD
  if(io_in->type == TMAP_SEQ_TYPE_SAM || io_in->type == TMAP_SEQ_TYPE_BAM) {
      // should be only one input file
      if(1 != io_in->n) {
          tmap_bug();
      }
      // get the current header
      if(NULL == io_in->seqios[0]) tmap_bug();
      if(NULL == io_in->seqios[0]->io.samio) tmap_bug();
      if(NULL == io_in->seqios[0]->io.samio->fp->header) tmap_bug();
      if(NULL == io_in->seqios[0]->io.samio->fp->header->header) {
          header = sam_header_parse2(io_in->seqios[0]->io.samio->fp->header->text);
      }
      else {
          header = io_in->seqios[0]->io.samio->fp->header->header; // wow, that's a lot of pointers
          if(NULL == header) tmap_bug();
          header = sam_header_clone(header); // clone the header
      }
      if(NULL == header) tmap_bug();
  }
  else {
      // empty header
      header = sam_header_init();
      // @HD - header line
      record = sam_header_record_init("HD"); // new header line
      if(0 == sam_header_record_add(record, "VN", "1.4")) tmap_bug(); // version number
      if(0 == sam_header_add_record(header, record)) tmap_bug(); // add the header line
      // nullify
      record = NULL;
  }

  // Get the TMAP program ID
  id = tmap_malloc(sizeof(char) * (1 + strlen(PACKAGE_NAME)), "id"); 
  strcpy(id, PACKAGE_NAME); // default
  for(i=j=0;NULL != (record_list = sam_header_get_record(header, "PG", "ID", id, &i)) && 0 < i;i=0) { // while the id is found
      char *ptr = NULL;
      // swap id and id_pp
      ptr = id_pp;
      id_pp = id;
      id = ptr;
      // create the new ID
      j++;
      id = tmap_realloc(id, sizeof(char) * (1 + (int)log10(j) + 1 + strlen(PACKAGE_NAME)), "id"); 
      if(sprintf(id, "%s.%d", PACKAGE_NAME, j) < 0) tmap_bug();
      free(record_list);
      record_list = NULL;
  }

  // @SQ
  if(NULL != refseq) {
      int32_t sq_exists = 0;
      sam_header_records_t *records = NULL;
      // NB: check to see if any SQ/SN records exist, if not, then ignore checking...
      records = sam_header_get_records(header, "SQ");
      if(NULL != records && 0 < records->n) sq_exists = 1;
      for(i=0;i<refseq->num_annos;i++) { // for each reference sequence
          char num[32];
          j = 0;
          if(1 == sq_exists) { // check to see if it already exists, if so ignore
              record_list = sam_header_get_record(header, "SQ", "SN", refseq->annos[i].name->s, &j);
          }
          if(0 == j) {
              record = sam_header_record_init("SQ"); // new reference sequence record
              if(0 == sam_header_record_add(record, "SN", refseq->annos[i].name->s)) tmap_bug(); // reference sequence name
              if(sprintf(num, "%u", (uint32_t)refseq->annos[i].len) < 0) tmap_bug(); // integer to string
              if(0 == sam_header_record_add(record, "LN", num)) tmap_bug(); // reference sequence length
              if(0 == sam_header_add_record(header, record)) tmap_bug(); // add the reference sequence record
          }
          else {
              free(record_list);
              record_list = NULL;
          }
      }
      record = NULL;
  }

  // @RG - read group
  if(0 < rg_sam_num) { // @RG specified on the command line
      // Check for SAM/BAM
      // TODO: this should be possible...
      if(io_in->type == TMAP_SEQ_TYPE_SAM || io_in->type == TMAP_SEQ_TYPE_BAM) {
          tmap_error("Cannot specify the read groups on the command line when using SAM/BAM as input."
                     "  Please embed in the SAM/BAM header instead.", Exit, OutOfRange);
      }
      record = NULL;
      // go through the command line arguments
      for(i=0;i<rg_sam_num;i++) {
          if(strlen(rg_sam[i]) < 4) tmap_error("Read group too small", Exit, OutOfRange);
          if(':' != rg_sam[i][2]) tmap_error("Read group improperly formatted (no colon)", Exit, OutOfRange);

          // check for id
          if('I' == rg_sam[i][0] && 'D' == rg_sam[i][1]) { // new read group
              if(NULL != record) { // add the record
                  tmap_seqs_io_init2_fs_and_add(io_in, header, record); // add @RG.KS and @RG.FO
              }
              record = sam_header_record_init("RG"); // new read group
          }
          // add the tag/value to the record
          if(NULL == record) {
              tmap_error("The read group ID must be specified first", Exit, OutOfRange);
          }
          tag[0]=rg_sam[i][0]; tag[1]=rg_sam[i][1]; // setup the tag
          if(0 == sam_header_record_add(record, tag, rg_sam[i]+3)) tmap_bug(); // add the tag/value
      }
      if(NULL != record) { // add the record
          tmap_seqs_io_init2_fs_and_add(io_in, header, record); // add @RG.KS and @RG.FO
      }
      // check that the # of read groups added was the same as the # of input files...
      sam_header_records_t *records = sam_header_get_records(header, "RG"); // get the header line
      if(records->n != io_in->n) tmap_error("The number of read groups did not match the number of input files", Exit, OutOfRange);
  }
  else if(io_in->type != TMAP_SEQ_TYPE_SAM && io_in->type != TMAP_SEQ_TYPE_BAM) { // dummy...
      for(i=0;i<io_in->n;i++) { // for each input file
          char buf[32];
          record = sam_header_record_init("RG"); // new read group
          if(1 == io_in->n) strcpy(buf, "NOID");
          else if(sprintf(buf, "NOID.%d", i+1) < 0) tmap_bug();
          if(0 == sam_header_record_add(record, "ID", buf)) tmap_bug(); // dummy ID
          if(0 == sam_header_record_add(record, "SM", "NOSM")) tmap_bug(); // dummy SM, for Picard validation
          if(0 == sam_header_record_add(record, "PG", id)) tmap_bug(); // dummy PG
          tmap_seqs_io_init2_fs_and_add(io_in, header, record); // add @RG.KS and @RG.FO
      }
  }
  else {
      // check that SM/PG are present
      sam_header_records_t *records = sam_header_get_records(header, "RG"); // get the header line
      for(i=0;i<records->n;i++) {
          record = records->records[i];
          if(NULL == sam_header_record_get(record, "ID")) tmap_error("Missing @RG.ID in the SAM/BAM Header", Exit, OutOfRange);
          if(NULL == sam_header_record_get(record, "SM")) {
              if(0 == sam_header_record_add(record, "SM", "NOSM")) tmap_bug(); // dummy SM, for Picard validation
          }
          if(NULL == sam_header_record_get(record, "PG")) {
              if(0 == sam_header_record_add(record, "PG", id)) tmap_bug(); // dummy PG
          }
      }
  }

  // @PG - program group
  // TODO: check for previous program group ID and set @PG.PP
  record = sam_header_record_init("PG"); // new program group
  if(0 == sam_header_record_add(record, "ID", id)) tmap_bug(); // @PG.ID
  if(0 == sam_header_record_add(record, "VN", PACKAGE_VERSION)) tmap_bug(); // @PG.VN
  // @PG.CL
  command_line = NULL;
  j = 1; // for the EOL
  command_line = tmap_realloc(command_line, sizeof(char) * j, "command_line");
  command_line[j-1] = '\0';
  for(i=0;i<argc;i++) {
      if(0 < i) j++;
      j += strlen(argv[i]);
      command_line = tmap_realloc(command_line, sizeof(char) * j, "command_line");
      if(0 < i) strcat(command_line, " ");
      strcat(command_line, argv[i]);
      command_line[j-1] = '\0';
  }
  if(0 == sam_header_record_add(record, "CL", command_line)) tmap_bug(); // @PG.CL
  if(NULL != id_pp) { // @PG.PP
      if(0 == sam_header_record_add(record, "PP", id_pp)) tmap_bug(); // @PG.CL
  }
  if(0 == sam_header_add_record(header, record)) tmap_bug(); // add the record
  free(command_line);

  // Check the new SAM Header
  if(0 == sam_header_check(header)) {
      tmap_error("SAM Header was not consistent", Exit, OutOfRange);
  }

  // Create a BAM Header from the SAM Header
  bam_header = bam_header_init(); // empty
  bam_header->header = header; // soft-copy the header
  bam_header = sam_header_to_bam_header(bam_header); // convert

  // free memory
  free(id);
  free(id_pp);

  return bam_header;
}

int
tmap_seqs_io_sff2sam_main(int argc, char *argv[])
{
  int c, help = 0;
  tmap_seqs_io_t *io_in = NULL;
  tmap_seqs_t *seqs = NULL;
  char **sam_rg = NULL;
  int32_t sam_rg_num = 0;
  int bidirectional = 0, sam_flowspace_tags = 0, remove_sff_clipping = 1;
  int out_type = 0;
  tmap_sam_io_t *io_out = NULL;
  bam_header_t *header = NULL; // BAM Header
  int32_t i;

  /*
  uint8_t *key_seq = NULL;
  int key_seq_len = 0;
  */

  while((c = getopt(argc, argv, "DGR:Yvh")) >= 0) {
      switch(c) {
        case 'D': bidirectional = 1; break;
        case 'G': remove_sff_clipping = 0; break;
        case 'R':
                  sam_rg = tmap_realloc(sam_rg, (1+sam_rg_num) * sizeof(char*), "sam_rg");
                  sam_rg[sam_rg_num] = tmap_strdup(optarg);
                  sam_rg_num++;
                  break;
        case 'Y': sam_flowspace_tags = 1; break;
        case 'v': tmap_progress_set_verbosity(1); break;
        case 'h': help = 1; break;
        default: return 1;
      }
  }
  if(1 != argc - optind || 1 == help) {
      tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [-R -Y -v -h] <in.sff>\n", PACKAGE, argv[0]);
      return 1; 
  }

  // input
  io_in = tmap_seqs_io_init(&argv[optind], 1, TMAP_SEQ_TYPE_SFF, TMAP_FILE_NO_COMPRESSION);

  // BAM Header
  header = tmap_seqs_io_to_bam_header(NULL, io_in, sam_rg, sam_rg_num, argc, argv);

  // open the output file
  switch(out_type) {
    case 0: // SAM
      io_out = tmap_sam_io_init2("-", "wh", header);
      break;
    case 1:
      io_out = tmap_sam_io_init2("-", "wb", header);
      break;
    case 2:
      io_out = tmap_sam_io_init2("-", "wbu", header);
      break;
    default:
      tmap_bug();
  }

  // destroy the BAM Header
  bam_header_destroy(header);
  header = NULL;

  seqs = tmap_seqs_init(TMAP_SEQ_TYPE_SFF);
  while(0 < tmap_seqs_io_read(io_in, seqs, io_out->fp->header->header)) {
      bam1_t *b = NULL;
      tmap_seq_t *seq = seqs->seqs[0];
      b = tmap_sam_convert_unmapped(seq, sam_flowspace_tags, bidirectional, NULL,
                                    0, 0, 0,
                                    0, 0, 0,
                                    "\tlq:i:%d\trq:i:%d\tla:i:%d\trq:i:%d",
                                    seq->data.sff->rheader->clip_qual_left,
                                    seq->data.sff->rheader->clip_qual_right,
                                    seq->data.sff->rheader->clip_adapter_left,
                                    seq->data.sff->rheader->clip_adapter_right);
      if(samwrite(io_out->fp, b) <= 0) {
          tmap_error("Error writing the SAM file", Exit, WriteFileError);
      }
      bam_destroy1(b); 
      tmap_seqs_destroy(seqs);
      seqs = tmap_seqs_init(TMAP_SEQ_TYPE_SFF);
  }
  tmap_seqs_destroy(seqs);

  // free memory
  tmap_seqs_io_destroy(io_in);
  tmap_sam_io_destroy(io_out);
  for(i=0;i<sam_rg_num;i++) {
      free(sam_rg[i]);
  }
  free(sam_rg);

  return 0;
}
