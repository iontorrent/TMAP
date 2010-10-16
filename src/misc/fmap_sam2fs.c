#include <stdlib.h>
#include <stdio.h>
#include <config.h>
#include <ctype.h>

#ifdef HAVE_SAMTOOLS
#include <sam.h>
#include <bam.h>
#endif

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_progress.h"
#include "../util/fmap_definitions.h"
#include "../io/fmap_file.h"
#include "../sw/fmap_fsw.h"
#include "fmap_sam2fs.h"

#ifdef HAVE_SAMTOOLS

// from bam.h
extern char *bam_nt16_rev_table;

// from fmap_fsw.c
extern int32_t fmap_fsw_sm_short[];

#define FMAP_SAM2FS_FLOW_MATCH '|'
#define FMAP_SAM2FS_FLOW_INS '+'
#define FMAP_SAM2FS_FLOW_DEL '-'
#define FMAP_SAM2FS_FLOW_SNP 'S'
#define FMAP_SAM2FS_FLOW_PAD ' '

static inline void
fmap_sam2fs_add_padding(char read_base, char *flow_order, int32_t *flow_j, 
                        uint8_t **read_flowgram, uint8_t **ref_flowgram, uint8_t **aln_flowgram, 
                        int32_t * flow_len, int32_t * flow_mem,
                        uint8_t aln_pad)
{
  int32_t k, padding = 0;

  while(flow_order[(*flow_j)] != read_base) {
      (*flow_j) = ((*flow_j)+1) & 0x3;
      padding++;
      if(4 <= padding) fmap_error("bug encountered", Exit, OutOfRange);
  }
  while((*flow_mem) <= (*flow_len) + padding + 1) { // padding 
      (*flow_mem) <<= 1;
      (*read_flowgram) = fmap_realloc((*read_flowgram), sizeof(uint8_t)*(*flow_mem), "read_flowgram");
      (*ref_flowgram) = fmap_realloc((*ref_flowgram), sizeof(uint8_t)*(*flow_mem), "ref_flowgram");
      (*aln_flowgram) = fmap_realloc((*aln_flowgram), sizeof(uint8_t)*(*flow_mem), "aln_flowgram");
  }
  for(k=0;k<padding;k++) {
      (*read_flowgram)[(*flow_len)+k] = (*ref_flowgram)[(*flow_len)+k] = 0;
      (*aln_flowgram)[(*flow_len)+k] = aln_pad;
  }
  (*flow_len) += padding;
}

static inline int32_t
fmap_sam2fs_is_DNA(char c)
{
  if('a' == c || 'A' == c
     || 'c' == c || 'C' == c
     || 'g' == c || 'G' == c
     || 't' == c || 'T' == c
     || 'n' == c || 'N' == c) {
      return 1;
  }
  else {
      return 0;
  }
}

static void 
fmap_sam2fs_aux(bam1_t *bam, char *flow_order, int32_t flow_score, int32_t flow_offset)
{
  int32_t i, j, k, l;

  uint8_t *md_data = NULL;
  char *md = NULL;
  uint32_t md_i = 0;

  uint32_t aln_len = 0;
  uint32_t *cigar = NULL;
  uint8_t *bam_seq = NULL;

  char *read_bases = NULL;
  int32_t read_bases_len = 0, read_bases_mem = 256;
  char *ref_bases = NULL;
  int32_t ref_bases_len = 0, ref_bases_mem = 64;
  uint8_t flow_order_tmp[4];

  int32_t soft_clip_start=0, soft_clip_end=0;

  int32_t flow_len;
  uint8_t *base_calls = NULL;
  uint16_t *flowgram = NULL;
  fmap_fsw_path_t *path = NULL;
  int32_t path_len;
  fmap_fsw_param_t param;
  int64_t score;

  // set the alignment parameters
  param.gap_open = 13*100;
  param.gap_ext = 2*100;
  param.gap_end = 2*100;
  param.matrix = fmap_fsw_sm_short;
  param.fscore = 100*flow_score;
  param.offset = flow_offset;
  param.row = 5;
  param.band_width = 50;

  if(BAM_FUNMAP & bam->core.flag) {
      return;
  }

  // get the MD tag
  if(NULL == (md_data = bam_aux_get(bam, "MD"))) {
      fmap_error("MD tag is missing", Exit, OutOfRange);
  }
  md = bam_aux2Z(md_data);

  // cigar
  cigar = bam1_cigar(bam);

  // get the alignment length
  for(i=aln_len=0;i<bam->core.n_cigar;i++) {
      aln_len += (cigar[i] >> 4);
  }
  if(0 == aln_len) {
      fmap_error("zero alignment length", Exit, OutOfRange);
  }

  // get soft clipping at the start/end of the sequence
  if((BAM_CSOFT_CLIP & cigar[0])) {
      soft_clip_start  = (cigar[0] >> 4);
  }
  if(1 < bam->core.n_cigar && (BAM_CSOFT_CLIP & cigar[bam->core.n_cigar-1])) {
      soft_clip_end = (cigar[bam->core.n_cigar-1] >> 4);
  }

  if(0 < soft_clip_start || 0 < soft_clip_end) {
      fmap_error("Soft clipping currently not supported", Exit, OutOfRange);
  }

  // get the read bases
  bam_seq = bam1_seq(bam);
  read_bases = fmap_calloc(1+bam->core.l_qseq, sizeof(char), "read_bases");
  for(i=0;i<bam->core.l_qseq;i++) {
      read_bases[i] = bam_nt16_rev_table[bam1_seqi(bam_seq, i)]; 
  }
  read_bases[i]='\0';
  read_bases_len = bam->core.l_qseq; 
  read_bases_mem = read_bases_len + 1;

  // get the reference bases
  ref_bases = fmap_calloc(ref_bases_mem, sizeof(char), "read_bases");

  // pre-process using the cigar array
  for(i=j=k=0;i<bam->core.n_cigar;i++) {
      int32_t op, op_len;
      op = (cigar[i] & 0xf);
      op_len = (cigar[i] >> 4);

      switch(op) {
        case BAM_CMATCH:
          // copy over these bases
          while(ref_bases_mem <= ref_bases_len + op_len + 1) {
              ref_bases_mem <<= 1;
              ref_bases = fmap_realloc(ref_bases, sizeof(char)*ref_bases_mem, "ref_bases");
          }
          for(l=0;l<op_len;l++) {
              ref_bases[k+l] = read_bases[j+l];
          }
          j += op_len; 
          k += op_len; 
          ref_bases_len += op_len;
          break;
        case BAM_CSOFT_CLIP:
        case BAM_CINS:
          // skip soft clipped and inserted bases 
          j += op_len;
          break;
        case BAM_CREF_SKIP:
        case BAM_CDEL:
        case BAM_CHARD_CLIP:
        case BAM_CPAD:
          // ignore
          break;
        default:
          fmap_error("unknown cigar operator", Exit, OutOfRange);
          break;
      }
  }

  // fill in with the MD array
  for(md_i=i=0;md_i<strlen(md);) {
      if('0' <= md[md_i] && md[md_i] <= '9') { // 0-9, matches
          l = atoi(md + md_i);
          i += l;
          while(md_i < strlen(md) && '0' <= md[md_i] && md[md_i] <= '9') { // 0-9
              md_i++; // skip over integers
          }
      }
      else if('^' == md[md_i]) { // deletion from the reference
          // how many bases are deleted?
          md_i++; // skip over '^'
          l=0;
          while(md_i+l < strlen(md) && 1 == fmap_sam2fs_is_DNA(md[md_i+l])) {
              l++;
          }
          // reallocate
          while(ref_bases_mem <= ref_bases_len + l + 1) { // more memory please
              ref_bases_mem <<= 1;
              ref_bases = fmap_realloc(ref_bases, sizeof(char)*ref_bases_mem, "ref_bases");
          }
          // shift up
          for(j=ref_bases_len-1;i<=j;j--) {
              ref_bases[j+l] = ref_bases[j];
          }
          // fill in
          for(j=0;j<l;j++) {
              ref_bases[i+j] = md[md_i+j];
          }
          md_i += l;
          i += l;
          ref_bases_len += l;
      }
      else if(1 == fmap_sam2fs_is_DNA(md[md_i])) { // SNP
          ref_bases[i] = md[md_i];
          i++;
          md_i++;
      }
      else {
          fmap_error("could not parse the MD tag", Exit, OutOfRange);
      }
  }
  ref_bases[i]='\0';
  //fprintf(stderr, "ref_bases_len=%d\nref_bases=%s\n", ref_bases_len, ref_bases);
  //fprintf(stderr, "read_bases_len=%d\nread_bases=%s\n", read_bases_len, read_bases);

  // DNA to integer
  for(i=0;i<read_bases_len;i++) {
      read_bases[i] = fmap_nt_char_to_int[(int)read_bases[i]];
  }
  for(i=0;i<ref_bases_len;i++) {
      ref_bases[i] = fmap_nt_char_to_int[(int)ref_bases[i]];
  }
  for(i=0;i<4;i++) {
      flow_order_tmp[i] = fmap_nt_char_to_int[(int)flow_order[i]];
  }

  flow_len = 0;
  for(i=j=0;i<read_bases_len;i++) {
      while(flow_order_tmp[j] != read_bases[i]) {
          flow_len++;
          j = (1+j) & 3;
      }
  }
  flow_len++; // the last flow

  base_calls = fmap_calloc(flow_len, sizeof(uint8_t), "base_calls");
  flowgram = fmap_calloc(flow_len, sizeof(uint16_t), "base_calls");

  for(i=j=0;i<flow_len;i++) {
      if(flow_order_tmp[i&3] == read_bases[j]) {
          k=j;
          while(j < read_bases_len 
                && read_bases[j] == read_bases[k]) {
              j++;
          }
          base_calls[i] = j-k;
          flowgram[i] = 100*(j-k);
      }
      else {
          base_calls[i] = 0;
          flowgram[i] = 0;
      }
  }

  path = fmap_calloc(FMAP_FSW_MAX_PATH_LENGTH(ref_bases_len, flow_len, param.offset), sizeof(fmap_fsw_path_t), "path"); 

  // re-align 
  score = fmap_fsw_global_core((uint8_t*)ref_bases, ref_bases_len,
                       flow_order_tmp, base_calls, flowgram, flow_len,
                       -1, 0,
                       &param, path, &path_len);

  // print
  fmap_fsw_print_aln(score, path, path_len, flow_order_tmp, (uint8_t*)ref_bases);

  // free memory
  free(path);
  free(base_calls);
  free(flowgram);
  free(read_bases);
  free(ref_bases);
}

static int
usage(char *flow_order, int32_t flow_score, int32_t flow_offset)
{
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Usage: %s sam2fs <in.sam/in.bam> [options]", PACKAGE);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (optional):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -f          the flow order [%s]\n", flow_order);
  fmap_file_fprintf(fmap_file_stderr, "         -F INT      the flow penalty [%d]\n", flow_score);
  fmap_file_fprintf(fmap_file_stderr, "         -o INT      search for homopolymer errors +- offset [%d]\n",
                    flow_offset);
  fmap_file_fprintf(fmap_file_stderr, "         -S          the input is a SAM file\n");
  fmap_file_fprintf(fmap_file_stderr, "         -v          print verbose progress information\n");
  fmap_file_fprintf(fmap_file_stderr, "         -h          print this message\n");
  fmap_file_fprintf(fmap_file_stderr, "\n");

  return 1;
}

static void
fmap_sam2fs_core(const char *fn_in, const char *sam_open_flags, char *flow_order,
                 int32_t flow_score, int32_t flow_offset)
{
  samfile_t *fp_in = NULL;
  bam1_t *b = NULL;

  fp_in = samopen(fn_in, sam_open_flags, 0);

  b = bam_init1();
  while(0 < samread(fp_in, b)) {
      // process
      fmap_sam2fs_aux(b, flow_order, flow_score, flow_offset);
      // destroy the bam
      bam_destroy1(b);
      // reinitialize
      b = bam_init1();
  }
  bam_destroy1(b);

  samclose(fp_in); 
}

int
fmap_sam2fs_main(int argc, char *argv[])
{
  int c;
  char sam_open_flags[16] = "rb";
  char *flow_order=NULL;
  int32_t flow_score, flow_offset;

  flow_order = fmap_strdup("TACG");
  flow_score = 26*100; // set this to score_match + gap_open + gap_ext
  flow_offset = 1;

  while((c = getopt(argc, argv, "f:F:o:Svh")) >= 0) {
      switch(c) {
        case 'f':
          strncpy(flow_order, optarg, 4); break;
        case 'F':
          flow_score = atoi(optarg); break;
        case 'o':
          flow_offset = atoi(optarg); break;
        case 'S':
          strcpy(sam_open_flags, "r"); break;
        case 'v':
          fmap_progress_set_verbosity(1); break;
        case 'h':
        default:
          return usage(flow_order, flow_score, flow_offset);
      }
  }


  if(argc != optind+1 || 1 == argc) {
      return usage(flow_order, flow_score, flow_offset);
  }
  else { // check command line options
      fmap_error_cmd_check_int(flow_score, 0, INT32_MAX, "-F");
      fmap_error_cmd_check_int(flow_offset, 0, INT32_MAX, "-o");
  }

  fmap_sam2fs_core(argv[optind], sam_open_flags, flow_order, flow_score, flow_offset);

  free(flow_order);

  fmap_progress_print2("terminating successfully");

  return 0;
}  
#endif
