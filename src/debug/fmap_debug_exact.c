#include <stdlib.h>
#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt_gen.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_bwt_match.h"
#include "../index/fmap_sa.h"
#include "../util/fmap_fq.h"
#include "../util/fmap_seq.h"
#include "../io/fmap_seq_io.h"
#include "fmap_debug_exact.h"

static void 
fmap_debug_exact_print_sam_unmapped(fmap_fq_t *seq)
{
  uint16_t flag = 0x0004;
  fmap_file_fprintf(fmap_file_stdout, "%s\t%u\t%s\t%u\t%u\t*\t*\t0\t0\t%s\t%s\n",
              seq->name->s, flag, "*",
              0, 0, seq->seq->s, seq->qual->s);
}

static int32_t 
fmap_debug_exact_print_sam(fmap_refseq_t *refseq, fmap_fq_t *seq, uint32_t pacpos, uint8_t strand)
{
  uint32_t pos, seqid;
  uint16_t flag = 0;

  if(0 < fmap_refseq_pac2real(refseq, pacpos, seq->seq->l, &seqid, &pos)) {
      if(1 == strand) { // reverse for the output
          flag |= 0x0010;
          fmap_string_reverse_compliment(seq->seq, seq->is_int);
          fmap_string_reverse(seq->qual);
      }
      fmap_file_fprintf(fmap_file_stdout, "%s\t%u\t%s\t%u\t%u\t%dM\t*\t0\t0\t%s\t%s\n",
              seq->name->s, flag, refseq->annos[seqid].name->s,
              pos, 255, (int)seq->seq->l, seq->seq->s, seq->qual->s);
      if(1 == strand) { // reverse back
          fmap_string_reverse_compliment(seq->seq, seq->is_int);
          fmap_string_reverse(seq->qual);
      }
      return 1;
  }
  return 0;
}

static void 
fmap_debug_exact_core_worker(fmap_refseq_t *refseq, fmap_bwt_t *bwt, fmap_sa_t *sa, fmap_fq_t *orig_seq)
{
  uint32_t i;
  uint32_t mapped = 0;
  fmap_fq_t *seq=NULL, *rseq=NULL;
  fmap_bwt_match_occ_t cur;

  seq = fmap_fq_clone(orig_seq);
  fmap_fq_to_int(seq);
  
  rseq = fmap_fq_clone(orig_seq);
  fmap_fq_to_int(rseq);
  fmap_string_reverse_compliment(rseq->seq, 1);

  if(0 != fmap_bwt_match_exact(bwt, seq->seq->l, (uint8_t*)seq->seq->s, &cur)) {
      for(i=cur.k;i<=cur.l;i++) {
          if(0 != fmap_debug_exact_print_sam(refseq, seq, fmap_sa_pac_pos(sa, bwt, i), 0)) {
              mapped = 1;
          }
      }
  }
  if(0 != fmap_bwt_match_exact(bwt, seq->seq->l, (uint8_t*)rseq->seq->s, &cur)) {
      for(i=cur.k;i<=cur.l;i++) {
          if(0 != fmap_debug_exact_print_sam(refseq, seq, fmap_sa_pac_pos(sa, bwt, i), 1)) {
              mapped = 1;
          }
      }
  }

  if(0 == mapped) {
      fmap_debug_exact_print_sam_unmapped(seq);
  }
  
  fmap_fq_destroy(seq);
  fmap_fq_destroy(rseq);
}

static void 
fmap_debug_exact_core(fmap_debug_exact_opt_t *opt)
{
  uint32_t i;
  fmap_refseq_t *refseq=NULL;
  fmap_bwt_t *bwt=NULL;
  fmap_sa_t *sa=NULL;
  fmap_file_t *fp_reads=NULL;
  fmap_seq_io_t *seqio=NULL;
  fmap_seq_t *seq=NULL;
  
  // SAM header
  refseq = fmap_refseq_read(opt->fn_fasta, 0);
  for(i=0;i<refseq->num_annos;i++) {
      fmap_file_fprintf(fmap_file_stdout, "@SQ\tSN:%s\tLN:%d\n",
              refseq->annos[i].name->s, (int)refseq->annos[i].len);
  }

  // TODO modify so we use forward search
  bwt = fmap_bwt_read(opt->fn_fasta, 1);
  sa = fmap_sa_read(opt->fn_fasta, 1);

  fp_reads = fmap_file_fopen(opt->fn_reads, "rb", FMAP_FILE_NO_COMPRESSION);
  seqio = fmap_seq_io_init(fp_reads, FMAP_SEQ_TYPE_FQ);
  seq = fmap_seq_init(FMAP_SEQ_TYPE_FQ);

  while(0 <= fmap_seq_io_read(seqio, seq)) {
      fmap_debug_exact_core_worker(refseq, bwt, sa, seq->data.fq);
  }

  fmap_file_fclose(fp_reads);
  fmap_refseq_destroy(refseq);
  fmap_bwt_destroy(bwt);
  fmap_sa_destroy(sa);
  fmap_seq_io_destroy(seqio);
  fmap_seq_destroy(seq);
}

static int 
usage(fmap_debug_exact_opt_t *opt)
{
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Usage: %s exact [options]", PACKAGE);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (required):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -f FILE     the FASTA reference file name\n");
  fmap_file_fprintf(fmap_file_stderr, "         -r FILE     the FASTQ reads file name\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (optional):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -h          print this message\n");
  fmap_file_fprintf(fmap_file_stderr, "\n");
  return 1;
}

int 
fmap_debug_exact(int argc, char *argv[])
{
  int c;
  fmap_debug_exact_opt_t opt;

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

  fmap_debug_exact_core(&opt);

  free(opt.fn_fasta);
  free(opt.fn_reads);

  return 0;
}
