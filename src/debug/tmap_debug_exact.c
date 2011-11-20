/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <unistd.h>
#include <config.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_sam_print.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt_gen.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_bwt_match.h"
#include "../index/tmap_sa.h"
#include "../seq/tmap_fq.h"
#include "../seq/tmap_seq.h"
#include "../seq/tmap_sam.h"
#include "../io/tmap_seq_io.h"
#include "tmap_debug_exact.h"

#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
static int32_t 
tmap_debug_exact_print_sam(tmap_refseq_t *refseq, tmap_fq_t *seq, tmp_bwt_int_t pacpos, uint8_t strand)
{
  tmp_bwt_int_t pos, seqid;
  uint16_t flag = 0;

  if(0 < tmap_refseq_pac2real(refseq, pacpos, seq->seq->l, &seqid, &pos)) {
      if(1 == strand) { // reverse for the output
          flag |= 0x0010;
          tmap_string_reverse_compliment(seq->seq, seq->is_int);
          tmap_string_reverse(seq->qual);
      }
      tmap_file_fprintf(tmap_file_stdout, "%s\t%u\t%s\t%u\t%u\t%dM\t*\t0\t0\t%s\t%s\n",
                        seq->name->s, flag, refseq->annos[seqid].name->s,
                        pos, 255, (int)seq->seq->l, seq->seq->s, seq->qual->s);
      if(1 == strand) { // reverse back
          tmap_string_reverse_compliment(seq->seq, seq->is_int);
          tmap_string_reverse(seq->qual);
      }
      return 1;
  }
  return 0;
}

static void 
tmap_debug_exact_core_worker(tmap_refseq_t *refseq, tmap_bwt_t *bwt, tmap_sa_t *sa, tmap_seq_t *orig_seq, int32_t n_only)
{
  uint32_t i, count = 0;
  uint32_t mapped = 0;
  tmp_bwt_int_t pacpos = 0;
  tmap_fq_t *seq=NULL, *rseq=NULL;
  tmap_bwt_match_occ_t cur;

  seq = tmap_fq_clone(orig_seq->data.fq);
  tmap_fq_to_int(seq);

  rseq = tmap_fq_clone(orig_seq->data.fq);
  tmap_fq_to_int(rseq);
  tmap_string_reverse_compliment(rseq->seq, 1);

  if(0 < tmap_bwt_match_exact(bwt, seq->seq->l, (uint8_t*)seq->seq->s, &cur)) {
      if(0 == n_only) {
          for(i=cur.k;i<=cur.l;i++) {
              pacpos = bwt->seq_len - tmap_sa_pac_pos(sa, bwt, i) - seq->seq->l + 1;
              if(0 != tmap_debug_exact_print_sam(refseq, orig_seq->data.fq, pacpos, 0)) {
                  mapped = 1;
              }
          }
      }
      else {
          count += cur.l - cur.k + 1;
      }
  }
  if(0 < tmap_bwt_match_exact(bwt, seq->seq->l, (uint8_t*)rseq->seq->s, &cur)) {
      if(0 == n_only) {
          for(i=cur.k;i<=cur.l;i++) {
              pacpos = bwt->seq_len - tmap_sa_pac_pos(sa, bwt, i) - seq->seq->l + 1;
              if(0 != tmap_debug_exact_print_sam(refseq, orig_seq->data.fq, pacpos, 1)) {
                  mapped = 1;
              }
          }
      }
      else {
          count += cur.l - cur.k + 1;
      }
  }

  if(0 == n_only) {
      if(0 == mapped) {
          tmap_sam_print_unmapped(tmap_file_stdout, orig_seq, 0, refseq, 0, 0, 0, 0, 0, 0);
      }
  }
  else {
      tmap_fq_to_char(seq);
      tmap_file_fprintf(tmap_file_stdout, "%s\t%s\t%u\n",
                        seq->name->s, seq->seq->s, count);
  }

  tmap_fq_destroy(seq);
  tmap_fq_destroy(rseq);
}

static void 
tmap_debug_exact_core(tmap_debug_exact_opt_t *opt)
{
  tmap_refseq_t *refseq=NULL;
  tmap_bwt_t *bwt=NULL;
  tmap_sa_t *sa=NULL;
  tmap_seq_io_t *seqio=NULL;
  tmap_seq_t *seq=NULL;
  
  seqio = tmap_seq_io_init(opt->fn_reads, 1, TMAP_SEQ_TYPE_FQ, TMAP_FILE_NO_COMPRESSION);
  seq = tmap_seq_init(TMAP_SEQ_TYPE_FQ);

  bwt = tmap_bwt_read(opt->fn_fasta, 1);
  if(0 == opt->n_only) {
      refseq = tmap_refseq_read(opt->fn_fasta, 0);
      sa = tmap_sa_read(opt->fn_fasta, 1);
      tmap_sam_print_header(tmap_file_stdout, refseq, seqio, NULL, NULL, NULL, 0, opt->argc, opt->argv);
  }

  while(0 <= tmap_seq_io_read(seqio, seq)) {
      tmap_debug_exact_core_worker(refseq, bwt, sa, seq, opt->n_only);
  }

  tmap_bwt_destroy(bwt);
  if(0 == opt->n_only) {
      tmap_refseq_destroy(refseq);
      tmap_sa_destroy(sa);
  }
  tmap_seq_io_destroy(seqio);
  tmap_seq_destroy(seq);
}

static int 
usage(tmap_debug_exact_opt_t *opt)
{
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Usage: %s exact [options]", PACKAGE);
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Options (required):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -f FILE     the FASTA reference file name\n");
  tmap_file_fprintf(tmap_file_stderr, "         -r FILE     the FASTQ reads file name\n");
  tmap_file_fprintf(tmap_file_stderr, "Options (optional):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -c          print only the number of hits per sequence\n");
  tmap_file_fprintf(tmap_file_stderr, "         -h          print this message\n");
  tmap_file_fprintf(tmap_file_stderr, "\n");
  return 1;
}

int 
tmap_debug_exact(int argc, char *argv[])
{
  int c;
  tmap_debug_exact_opt_t opt;

  opt.fn_fasta = NULL;
  opt.fn_reads = NULL;
  opt.n_only = 0;

  while((c = getopt(argc, argv, "f:r:ch")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = tmap_strdup(optarg); break;
        case 'r':
          opt.fn_reads = tmap_strdup(optarg); break;
        case 'c':
          opt.n_only = 1; break;
        case 'h':
        default:
          return usage(&opt);
      }
  }

  if(argc != optind || 1 == argc) {
      return usage(&opt);
  }
  if(NULL == opt.fn_fasta) {
      tmap_error("required option -f", Exit, CommandLineArgument);
  }
  if(NULL == opt.fn_reads) {
      tmap_error("required option -r", Exit, CommandLineArgument);
  }

  // Note: 'tmap_file_stdout' should not have been previously modified
  tmap_file_stdout = tmap_file_fdopen(fileno(stdout), "wb", TMAP_FILE_NO_COMPRESSION);

  tmap_debug_exact_core(&opt);

  // close the output
  tmap_file_fclose(tmap_file_stdout);

  free(opt.fn_fasta);
  free(opt.fn_reads);

  return 0;
}
#endif
