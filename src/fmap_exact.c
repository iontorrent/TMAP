#include <stdlib.h>
#include "fmap_error.h"
#include "fmap_alloc.h"
#include "fmap_refseq.h"
#include "fmap_bwt_gen.h"
#include "fmap_bwt.h"
#include "fmap_sa.h"
#include "fmap_seq.h"
#include "fmap_definitions.h"
#include "fmap_exact.h"

static void fmap_exact_print_sam_unmapped(fmap_seq_t *seq)
{
  uint16_t flag = 0x0004;
  fprintf(stdout, "%s\t%u\t%s\t%u\t%u\t*\t*\t0\t0\t%s\t%s\n",
              seq->name.s, flag, "*",
              0, 0, seq->seq.s, seq->qual.s);
}

static int32_t fmap_exact_print_sam(fmap_refseq_t *refseq, fmap_seq_t *seq, uint32_t pacpos, uint8_t strand)
{
  uint32_t pos, seqid;
  uint16_t flag = 0;

  if(0 <= fmap_refseq_pac2real(refseq, pacpos, seq->seq.l, &seqid, &pos)) {
      if(1 == strand) { // reverse for the outptu
          flag |= 0x0010;
          fmap_seq_reverse(seq, 1);
      }
      fprintf(stdout, "%s\t%u\t%s\t%u\t%u\t%dM\t*\t0\t0\t%s\t%s\n",
              seq->name.s, flag, refseq->annos[seqid].name,
              pos, 255, (int)seq->seq.l, seq->seq.s, seq->qual.s);
      if(1 == strand) { // reverse back
          fmap_seq_reverse(seq, 1);
      }
      return 1;
  }
  return 0;
}

static void fmap_exact_core_worker(fmap_refseq_t *refseq, fmap_bwt_t *bwt, fmap_sa_t *sa, fmap_seq_t *seq)
{
  uint32_t i;
  uint32_t sa_begin, sa_end;
  uint8_t *seq_int = NULL, *seq_int_rc = NULL;
  uint32_t mapped = 0;

  seq_int = fmap_malloc(sizeof(uint8_t)*seq->seq.l, "seq_int");
  seq_int_rc = fmap_malloc(sizeof(uint8_t)*seq->seq.l, "seq_int");
  for(i=0;i<seq->seq.l;i++) {
      seq_int[i] = nt_char_to_int[(int)seq->seq.s[i]];
      seq_int_rc[seq->seq.l-i-1] = (4 <= seq_int[i]) ? seq_int[i] : 3 - seq_int[i];
  }

  if(0 != bwt_match_exact(bwt, seq->seq.l, seq_int, &sa_begin, &sa_end)) {
      for(i=sa_begin;i<=sa_end;i++) {
          if(0 != fmap_exact_print_sam(refseq, seq, fmap_sa_pac_pos(sa, bwt, i), 0)) {
              mapped = 1;
          }
      }
  }
  if(0 != bwt_match_exact(bwt, seq->seq.l, seq_int_rc, &sa_begin, &sa_end)) {
      for(i=sa_begin;i<=sa_end;i++) {
          if(0 != fmap_exact_print_sam(refseq, seq, fmap_sa_pac_pos(sa, bwt, i), 1)) {
              mapped = 1;
          }
      }
  }

  if(0 == mapped) {
      fmap_exact_print_sam_unmapped(seq);
  }
  
  free(seq_int);
}

static void fmap_exact_core(fmap_exact_opt_t *opt)
{
  uint32_t i;
  fmap_refseq_t *refseq=NULL;
  fmap_bwt_t *bwt=NULL;
  fmap_sa_t *sa=NULL;
  fmap_file_t *fp_reads=NULL;
  fmap_seq_t *seq=NULL;
  
  // SAM header
  refseq = fmap_refseq_read(opt->fn_fasta);
  for(i=0;i<refseq->num_annos;i++) {
      fprintf(stdout, "@SQ\tSN:%s\tLN:%d\n",
              refseq->annos[i].name, (int)refseq->annos[i].len);
  }

  bwt = fmap_bwt_read(opt->fn_fasta);
  sa = fmap_sa_read(opt->fn_fasta);

  fp_reads = fmap_file_fopen(opt->fn_reads, "rb", FMAP_FILE_NO_COMPRESSION);
  seq = fmap_seq_init(fp_reads);

  while(0 <= fmap_seq_read(seq)) {
      fmap_exact_core_worker(refseq, bwt, sa, seq);
  }

  fmap_file_fclose(fp_reads);
  fmap_refseq_destroy(refseq);
  fmap_bwt_destroy(bwt);
  fmap_sa_destroy(sa);
}

static int usage(fmap_exact_opt_t *opt)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage: %s exact [optionsn", PACKAGE);
  fprintf(stderr, "\n");
  fprintf(stderr, "Optiosn (required):\n");
  fprintf(stderr, "         -f FILE     the FASTA reference file name\n");
  fprintf(stderr, "         -r FILE     the FASTQ reads file name\n");
  fprintf(stderr, "Options (optional):\n");
  fprintf(stderr, "         -h          print this message\n");
  fprintf(stderr, "\n");
  return 1;
}

int fmap_exact(int argc, char *argv[])
{
  int c;
  fmap_exact_opt_t opt;

  opt.fn_fasta = NULL;

  while((c = getopt(argc, argv, "f:r:h")) >= 0) {
      switch(c) {
        case 'f':
          opt.fn_fasta = fmap_strdup(optarg); break;
        case 'r':
          opt.fn_reads = fmap_strdup(optarg); break;
        case 'h':
        default:
          return usage(&opt);
      }
  }

  if(argc != optind || 1 == argc) {
      return usage(&opt);
  }
  if(NULL == opt.fn_fasta) {
      fmap_error("required option -f", Exit, CommandLineArgument);
  }
  if(NULL == opt.fn_reads) {
      fmap_error("required option -r", Exit, CommandLineArgument);
  }

  fmap_exact_core(&opt);

  free(opt.fn_fasta);
  free(opt.fn_reads);

  return 0;
}
