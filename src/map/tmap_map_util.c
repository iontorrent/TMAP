/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "../util/tmap_alloc.h"
#include "../util/tmap_error.h"
#include "../util/tmap_sam_print.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sort.h"
#include "../util/tmap_definitions.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_sa.h"
#include "../sw/tmap_sw.h"
#include "../sw/tmap_fsw.h"
#include "tmap_map_util.h"

#define tmap_map_util_reverse_query(_query, _ql, _i) \
  for(_i=0;_i<(_ql>>1);_i++) { \
      uint8_t _tmp = _query[_i]; \
      _query[_i] = _query[_ql-_i-1]; \
      _query[_ql-_i-1] = _tmp; \
  }

// sort by strand, min-seqid, min-position
#define __tmap_map_sam_sort_coord_lt(a, b) (  ((a).strand < (b).strand) \
                                            || ( (a).strand == (b).strand && (a).seqid < (b).seqid) \
                                            || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos < (b).pos ) \
                                            ? 1 : 0 )

// sort by strand, min-seqid, min-position, max score
#define __tmap_map_sam_sort_coord_score_lt(a, b) (  ((a).strand < (b).strand) \
                                            || ( (a).strand == (b).strand && (a).seqid < (b).seqid) \
                                            || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos < (b).pos ) \
                                            || ( (a).strand == (b).strand && (a).seqid == (b).seqid && (a).pos == (b).pos && (a).score > (b).score) \
                                            ? 1 : 0 )

// sort by max-score
#define __tmap_map_sam_sort_score_lt(a, b) ((a).score > (b).score)

TMAP_SORT_INIT(tmap_map_sam_sort_coord, tmap_map_sam_t, __tmap_map_sam_sort_coord_lt)
TMAP_SORT_INIT(tmap_map_sam_sort_coord_score, tmap_map_sam_t, __tmap_map_sam_sort_coord_score_lt)
TMAP_SORT_INIT(tmap_map_sam_sort_score, tmap_map_sam_t, __tmap_map_sam_sort_score_lt)

tmap_map_opt_t *
tmap_map_opt_init(int32_t algo_id)
{
  int32_t i;
  tmap_map_opt_t *opt = NULL;

  opt = tmap_calloc(1, sizeof(tmap_map_opt_t), "opt");

  // internal data
  opt->algo_id = algo_id;
  opt->algo_stage = -1;
  opt->argv = NULL;
  opt->argc = -1;

  // global options 
  opt->fn_fasta = opt->fn_reads = NULL;
  opt->reads_format = TMAP_READS_FORMAT_UNKNOWN;
  opt->score_match = TMAP_MAP_UTIL_SCORE_MATCH;
  opt->pen_mm = TMAP_MAP_UTIL_PEN_MM;
  opt->pen_gapo = TMAP_MAP_UTIL_PEN_GAPO;
  opt->pen_gape = TMAP_MAP_UTIL_PEN_GAPE;
  opt->fscore = TMAP_MAP_UTIL_FSCORE;
  opt->flow_order = NULL;
  opt->flow_order_use_sff = 0;
  opt->bw = 50; 
  opt->softclip_type = TMAP_MAP_UTIL_SOFT_CLIP_RIGHT;
  opt->dup_window = 128;
  opt->max_seed_band = 15;
  opt->score_thr = 8;
  opt->reads_queue_size = 262144;
  opt->num_threads = 1;
  opt->aln_output_mode = TMAP_MAP_UTIL_ALN_MODE_RAND_BEST;
  opt->sam_rg = NULL;
  opt->sam_sff_tags = 0;
  opt->input_compr = TMAP_FILE_NO_COMPRESSION;
  opt->output_compr = TMAP_FILE_NO_COMPRESSION;
  opt->shm_key = 0;
  opt->min_seq_len = -1;
  opt->max_seq_len = -1;

  switch(algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      // map1
      opt->seed_length = 32;
      opt->seed_max_diff = 2;
      opt->seed2_length = 64;
      opt->max_diff = -1; opt->max_diff_fnr = 0.04;
      opt->max_mm = 3; opt->max_mm_frac = -1.;
      opt->max_err_rate = 0.02;
      opt->max_gapo = 1; opt->max_gapo_frac = -1.;
      opt->max_gape = 6; opt->max_gape_frac = -1.;
      opt->max_cals_del = 10;
      opt->indel_ends_bound = 5;
      opt->max_best_cals = 32;
      opt->max_entries = 2000000;
      break;
    case TMAP_MAP_ALGO_MAP2:
      // map2
      opt->yita = 5.5f;
      //opt->mask_level = 0.50; 
      opt->length_coef = 5.5f;
      opt->max_seed_intv = 3; 
      opt->z_best = 1; 
      opt->seeds_rev = 5;
      break;
    case TMAP_MAP_ALGO_MAP3:
      // map3
      opt->seed_length = -1;
      opt->seed_length_set = 0;
      opt->max_seed_hits = 12;
      opt->hp_diff = 0;
      break;
    case TMAP_MAP_ALGO_MAPALL:
      // mapall
      opt->aln_output_mode_ind = 0;
      for(i=0;i<2;i++) {
          opt->algos[i] = 0;
          opt->opt_map1[i] = tmap_map_opt_init(TMAP_MAP_ALGO_MAP1);
          opt->opt_map2[i] = tmap_map_opt_init(TMAP_MAP_ALGO_MAP2);
          opt->opt_map3[i] = tmap_map_opt_init(TMAP_MAP_ALGO_MAP3);
      }
      break;
    default:
      break;
  }

  return opt;
}

void
tmap_map_opt_destroy(tmap_map_opt_t *opt)
{
  int32_t i;

  free(opt->fn_fasta);
  free(opt->fn_reads); 
  free(opt->sam_rg);
  free(opt->flow_order);

  switch(opt->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
    case TMAP_MAP_ALGO_MAP2:
    case TMAP_MAP_ALGO_MAP3:
      break;
    case TMAP_MAP_ALGO_MAPALL:
      // mapall
      for(i=0;i<2;i++) {
          opt->algos[i] = 0;
          tmap_map_opt_destroy(opt->opt_map1[i]);
          tmap_map_opt_destroy(opt->opt_map2[i]);
          tmap_map_opt_destroy(opt->opt_map3[i]);
      }
      break;
    default:
      break;
  }

  free(opt);
}

#define __tmap_map_print_compression(_type) switch(_type) { \
  case TMAP_FILE_NO_COMPRESSION: \
                                 tmap_file_fprintf(tmap_file_stderr, " [none]\n"); \
    break; \
  case TMAP_FILE_GZ_COMPRESSION: \
                                 tmap_file_fprintf(tmap_file_stderr, " [gz]\n"); \
    break; \
  case TMAP_FILE_BZ2_COMPRESSION: \
                                  tmap_file_fprintf(tmap_file_stderr, " [bz2]\n"); \
    break; \
  default: \
           tmap_file_fprintf(tmap_file_stderr, " [?]\n"); \
    break; \
}

#define __tmap_map_opt_usage_map1(_opt, _stage) do { \
    if(_stage < 0) tmap_file_fprintf(tmap_file_stderr, "%s options (optional):\n", tmap_algo_id_to_name(TMAP_MAP_ALGO_MAP1)); \
    else tmap_file_fprintf(tmap_file_stderr, "%s stage %d options (optional):\n", tmap_algo_id_to_name(TMAP_MAP_ALGO_MAP1), _stage); \
    tmap_file_fprintf(tmap_file_stderr, "         -l INT      the k-mer length to seed CALs (-1 to disable) [%d]\n", (_opt)->seed_length); \
    tmap_file_fprintf(tmap_file_stderr, "         -s INT      maximum number of edits in the seed [%d]\n", (_opt)->seed_max_diff); \
    tmap_file_fprintf(tmap_file_stderr, "         -L INT      the secondary seed length (-1 to disable) [%d]\n", (_opt)->seed2_length); \
    tmap_file_fprintf(tmap_file_stderr, "         -p NUM      maximum number of edits or false-negative probability assuming the maximum error rate "); \
    if((_opt)->max_diff < 0) tmap_file_fprintf(tmap_file_stderr, "[number: %d]\n", (_opt)->max_diff); \
    else tmap_file_fprintf(tmap_file_stderr, "[probability: %d]\n", (_opt)->max_diff_fnr); \
    tmap_file_fprintf(tmap_file_stderr, "         -P NUM      the assumed per-base maximum error rate [%lf]\n", (_opt)->max_err_rate); \
    tmap_file_fprintf(tmap_file_stderr, "         -m NUM      maximum number of or (read length) fraction of mismatches"); \
    if((_opt)->max_mm < 0) tmap_file_fprintf(tmap_file_stderr, " [fraction: %lf]\n", (_opt)->max_mm_frac); \
    else tmap_file_fprintf(tmap_file_stderr, " [number: %d]\n", (_opt)->max_mm); \
    tmap_file_fprintf(tmap_file_stderr, "         -o NUM      maximum number of or (read length) fraction of indel starts"); \
    if((_opt)->max_gapo < 0) tmap_file_fprintf(tmap_file_stderr, " [fraction: %lf]\n", (_opt)->max_gapo_frac); \
    else tmap_file_fprintf(tmap_file_stderr, " [number: %d]\n", (_opt)->max_gapo); \
    tmap_file_fprintf(tmap_file_stderr, "         -e NUM      maximum number of or (read length) fraction of indel extensions"); \
    if((_opt)->max_gape < 0) tmap_file_fprintf(tmap_file_stderr, " [fraction: %lf]\n", (_opt)->max_gape_frac); \
    else tmap_file_fprintf(tmap_file_stderr, " [number: %d]\n", (_opt)->max_gape); \
    tmap_file_fprintf(tmap_file_stderr, "         -d INT      the maximum number of CALs to extend a deletion [%d]\n", (_opt)->max_cals_del); \
    tmap_file_fprintf(tmap_file_stderr, "         -i INT      indels are not allowed within INT number of bps from the end of the read [%d]\n", (_opt)->indel_ends_bound); \
    tmap_file_fprintf(tmap_file_stderr, "         -b INT      stop searching when INT optimal CALs have been found [%d]\n", (_opt)->max_best_cals); \
    tmap_file_fprintf(tmap_file_stderr, "         -Q INT      maximum number of alignment nodes [%d]\n", (_opt)->max_entries); \
    tmap_file_fprintf(tmap_file_stderr, "         -u INT      the minimum sequence length to examine [%d]\n", (_opt)->min_seq_len); \
    tmap_file_fprintf(tmap_file_stderr, "         -U INT      the maximum sequence length to examine [%d]\n", (_opt)->max_seq_len); \
    tmap_file_fprintf(tmap_file_stderr, "\n"); \
    } while(0)

#define __tmap_map_opt_usage_map2(_opt, _stage) do { \
    if(_stage < 0) tmap_file_fprintf(tmap_file_stderr, "%s options (optional):\n", tmap_algo_id_to_name(TMAP_MAP_ALGO_MAP2)); \
    else tmap_file_fprintf(tmap_file_stderr, "%s stage %d options (optional):\n", tmap_algo_id_to_name(TMAP_MAP_ALGO_MAP2), _stage); \
    tmap_file_fprintf(tmap_file_stderr, "         -c FLOAT    coefficient of length-threshold adjustment [%.1lf]\n", (_opt)->length_coef); \
    tmap_file_fprintf(tmap_file_stderr, "         -S INT      maximum seeding interval size [%d]\n", (_opt)->max_seed_intv); \
    tmap_file_fprintf(tmap_file_stderr, "         -b INT      Z-best [%d]\n", (_opt)->z_best); \
    tmap_file_fprintf(tmap_file_stderr, "         -N INT      # seeds to trigger reverse alignment [%d]\n", (_opt)->seeds_rev); \
    tmap_file_fprintf(tmap_file_stderr, "         -u INT      the minimum sequence length to examine [%d]\n", (_opt)->min_seq_len); \
    tmap_file_fprintf(tmap_file_stderr, "         -U INT      the maximum sequence length to examine [%d]\n", (_opt)->max_seq_len); \
    tmap_file_fprintf(tmap_file_stderr, "\n"); \
    } while(0)

#define __tmap_map_opt_usage_map3(_opt, _stage) do { \
    if(_stage < 0) tmap_file_fprintf(tmap_file_stderr, "%s options (optional):\n", tmap_algo_id_to_name(TMAP_MAP_ALGO_MAP3)); \
    else tmap_file_fprintf(tmap_file_stderr, "%s stage %d options (optional):\n", tmap_algo_id_to_name(TMAP_MAP_ALGO_MAP3), _stage); \
    tmap_file_fprintf(tmap_file_stderr, "         -l INT      the k-mer length to seed CALs (-1 tunes to the genome size) [%d]\n", (_opt)->seed_length); \
    tmap_file_fprintf(tmap_file_stderr, "         -S INT      the maximum number of hits returned by a seed [%d]\n", (_opt)->max_seed_hits); \
    tmap_file_fprintf(tmap_file_stderr, "         -H INT      single homopolymer error difference for enumeration [%d]\n", (_opt)->hp_diff); \
    tmap_file_fprintf(tmap_file_stderr, "         -u INT      the minimum sequence length to examine [%d]\n", (_opt)->min_seq_len); \
    tmap_file_fprintf(tmap_file_stderr, "         -U INT      the maximum sequence length to examine [%d]\n", (_opt)->max_seq_len); \
    tmap_file_fprintf(tmap_file_stderr, "\n"); \
    } while(0)

int
tmap_map_opt_usage(tmap_map_opt_t *opt)
{
  char *reads_format = tmap_get_reads_file_format_string(opt->reads_format);
  int32_t i;

  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Usage: %s %s [options]\n", tmap_algo_id_to_name(opt->algo_id), PACKAGE);
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "global options (required):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -f FILE     the FASTA reference file name [%s]\n", opt->fn_fasta);
  tmap_file_fprintf(tmap_file_stderr, "         -r FILE     the reads file name [%s]\n", (NULL == opt->fn_reads) ? "stdin" : opt->fn_reads);
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "global options (optional):\n");
#ifdef HAVE_SAMTOOLS
  tmap_file_fprintf(tmap_file_stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff|sam|bam) [%s]\n", reads_format);
#else
  tmap_file_fprintf(tmap_file_stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff) [%s]\n", reads_format);
#endif
  tmap_file_fprintf(tmap_file_stderr, "         -A INT      score for a match [%d]\n", opt->score_match);
  tmap_file_fprintf(tmap_file_stderr, "         -M INT      the mismatch penalty [%d]\n", opt->pen_mm);
  tmap_file_fprintf(tmap_file_stderr, "         -O INT      the indel start penalty [%d]\n", opt->pen_gapo);
  tmap_file_fprintf(tmap_file_stderr, "         -E INT      the indel extend penalty [%d]\n", opt->pen_gape);
  tmap_file_fprintf(tmap_file_stderr, "         -X INT      the flow score penalty [%d]\n", opt->fscore);
  tmap_file_fprintf(tmap_file_stderr, "         -x STRING   the flow order ([ACGT]{4+} or \"sff\") [%s]\n",
                    (NULL == opt->flow_order) ? "not using" : opt->flow_order);
  tmap_file_fprintf(tmap_file_stderr, "         -w INT      the band width [%d]\n", opt->bw);
  tmap_file_fprintf(tmap_file_stderr, "         -g          the soft-clipping type [%d]\n", opt->softclip_type);
  tmap_file_fprintf(tmap_file_stderr, "                             0 - allow on the right and left portions of the read\n");
  tmap_file_fprintf(tmap_file_stderr, "                             1 - allow on the left portion of the read\n");
  tmap_file_fprintf(tmap_file_stderr, "                             2 - allow on the right portion of the read\n");
  tmap_file_fprintf(tmap_file_stderr, "                             3 - do not allow soft-clipping\n");
  tmap_file_fprintf(tmap_file_stderr, "         -W INT      remove duplicate alignments within this bp window (-1 to disable) [%d]\n",
                    opt->dup_window);
  tmap_file_fprintf(tmap_file_stderr, "         -B INT      the window of bases in which to group seeds [%d]\n", opt->max_seed_band); 
  tmap_file_fprintf(tmap_file_stderr, "         -T INT      score threshold divided by the match score [%d]\n", opt->score_thr);
  tmap_file_fprintf(tmap_file_stderr, "         -q INT      the queue size for the reads (-1 disables) [%d]\n", opt->reads_queue_size);
  tmap_file_fprintf(tmap_file_stderr, "         -n INT      the number of threads [%d]\n", opt->num_threads);
  tmap_file_fprintf(tmap_file_stderr, "         -a INT      output filter [%d]\n", opt->aln_output_mode);
  tmap_file_fprintf(tmap_file_stderr, "                             0 - unique best hits\n");
  tmap_file_fprintf(tmap_file_stderr, "                             1 - random best hit\n");
  tmap_file_fprintf(tmap_file_stderr, "                             2 - all best hits\n");
  tmap_file_fprintf(tmap_file_stderr, "                             3 - all alignments\n");
  tmap_file_fprintf(tmap_file_stderr, "         -R STRING   the RG tags to add to the SAM header [%s]\n", opt->sam_rg);
  tmap_file_fprintf(tmap_file_stderr, "         -Y          include SFF specific SAM tags [%s]\n",
                    (1 == opt->sam_sff_tags) ? "true" : "false");
  tmap_file_fprintf(tmap_file_stderr, "         -z/-j       the input is gz/bz2 compressed (gzip/bzip2)");
  __tmap_map_print_compression(opt->input_compr);
  tmap_file_fprintf(tmap_file_stderr, "         -Z/-J       the output is gz/bz2 compressed (gzip/bzip2)");
  __tmap_map_print_compression(opt->output_compr);
  tmap_file_fprintf(tmap_file_stderr, "         -k INT      use shared memory with the following key [%d]\n", opt->shm_key);
  tmap_file_fprintf(tmap_file_stderr, "         -v          print verbose progress information\n");
  tmap_file_fprintf(tmap_file_stderr, "         -h          print this message\n");
  tmap_file_fprintf(tmap_file_stderr, "\n");

  switch(opt->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      __tmap_map_opt_usage_map1(opt, -1);
      break;
    case TMAP_MAP_ALGO_MAP2:
      __tmap_map_opt_usage_map2(opt, -1);
      break;
    case TMAP_MAP_ALGO_MAP3:
      __tmap_map_opt_usage_map3(opt, -1);
      break;
    case TMAP_MAP_ALGO_MAPALL:
      for(i=0;i<2;i++) {
          __tmap_map_opt_usage_map1(opt->opt_map1[i], i+1);
          __tmap_map_opt_usage_map2(opt->opt_map2[i], i+1);
          __tmap_map_opt_usage_map3(opt->opt_map3[i], i+1);
      }
      // mapall
      tmap_file_fprintf(tmap_file_stderr, "%s options (optional):\n", tmap_algo_id_to_name(opt->algo_id));
      tmap_file_fprintf(tmap_file_stderr, "         -I          apply the output filter (-a) and duplicate removal (-W) for each algorithm separately [%s]\n",
                        (1 == opt->aln_output_mode_ind) ? "true" : "false");
      tmap_file_fprintf(tmap_file_stderr, "\n");
      break;
    default:
      break;
  }

  // free
  free(reads_format);

  return 1;
}

int32_t
tmap_map_opt_parse(int argc, char *argv[], tmap_map_opt_t *opt)
{
  int c;
  char *getopt_format = NULL;

  opt->argc = argc; opt->argv = argv;
  switch(opt->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      getopt_format = tmap_strdup("f:r:F:A:M:O:E:X:x:w:g:W:B:T:q:n:a:R:YjzJZk:vhl:s:L:p:P:m:o:e:d:i:b:Q:u:U:");
      break;
    case TMAP_MAP_ALGO_MAP2:
      getopt_format = tmap_strdup("f:r:F:A:M:O:E:X:x:w:g:W:B:T:q:n:a:R:YjzJZk:vhc:S:b:N:u:U:");
      break;
    case TMAP_MAP_ALGO_MAP3:
      getopt_format = tmap_strdup("f:r:F:A:M:O:E:X:x:w:g:W:B:T:q:n:a:R:YjzJZk:vhl:S:H:u:U:");
      break;
    case TMAP_MAP_ALGO_MAPALL:
      getopt_format = tmap_strdup("f:r:F:A:M:O:E:X:x:w:g:W:B:T:q:n:a:R:YjzJZk:vhIu:U:");
      break;
    default:
      break;
  }

  // global options
  // Note: possible memory leaks if the same option (besides -R) are specified twice
  while((c = getopt(argc, argv, getopt_format)) >= 0) {
      switch(c) { 
        case 'f': 
          opt->fn_fasta = tmap_strdup(optarg); break;
        case 'r':
          opt->fn_reads = tmap_strdup(optarg); 
          tmap_get_reads_file_format_from_fn_int(opt->fn_reads, &opt->reads_format, &opt->input_compr);
          break;
        case 'F':
          opt->reads_format = tmap_get_reads_file_format_int(optarg); break;
        case 'A':
          opt->score_match = atoi(optarg); break;
        case 'M':
          opt->pen_mm = atoi(optarg); break;
        case 'O':
          opt->pen_gapo = atoi(optarg); break;
        case 'E':
          opt->pen_gape = atoi(optarg); break;
        case 'X':
          opt->fscore = atoi(optarg); break;
        case 'x':
          opt->flow_order = tmap_strdup(optarg); 
          if(0 == strcmp("sff", opt->flow_order) || 0 == strcmp("SFF", opt->flow_order)) {
              opt->flow_order_use_sff = 1;
          }
          else {
              opt->flow_order_use_sff = 0;
          }
          break;
        case 'w':
          opt->bw = atoi(optarg); break;
        case 'g':
          opt->softclip_type = atoi(optarg); break;
        case 'W':
          opt->dup_window = atoi(optarg); break;
        case 'B':
          opt->max_seed_band = atoi(optarg); break;
        case 'T':
          opt->score_thr = atoi(optarg); break;
        case 'q':
          opt->reads_queue_size = atoi(optarg); break;
        case 'n':
          opt->num_threads = atoi(optarg); break;
        case 'a':
          opt->aln_output_mode = atoi(optarg); break;
        case 'R':
          if(NULL == opt->sam_rg) {
              // add fiv for the string "@RG\t" and null terminator
              opt->sam_rg = tmap_realloc(opt->sam_rg, sizeof(char) * (5 + strlen(optarg)), "opt->sam_rg");
              strcpy(opt->sam_rg, "@RG\t");
              strcat(opt->sam_rg, optarg);
          }
          else {
              // add two for the tab separator and null terminator
              opt->sam_rg = tmap_realloc(opt->sam_rg, sizeof(char) * (2 + strlen(optarg) + strlen(opt->sam_rg)), "opt->sam_rg");
              if(0 < strlen(optarg) && '\t' != optarg[0]) strcat(opt->sam_rg, "\t"); // add a tab separator
              strcat(opt->sam_rg, optarg);
          }
          // remove trailing white spaces
          tmap_chomp(opt->sam_rg);
          break;
        case 'Y':
          opt->sam_sff_tags = 1; break;
        case 'j':
          opt->input_compr = TMAP_FILE_BZ2_COMPRESSION;
          tmap_get_reads_file_format_from_fn_int(opt->fn_reads, &opt->reads_format, &opt->input_compr);
          break;
        case 'z':
          opt->input_compr = TMAP_FILE_GZ_COMPRESSION;
          tmap_get_reads_file_format_from_fn_int(opt->fn_reads, &opt->reads_format, &opt->input_compr);
          break;
        case 'J':
          opt->output_compr = TMAP_FILE_BZ2_COMPRESSION; break;
        case 'Z':
          opt->output_compr = TMAP_FILE_GZ_COMPRESSION; break;
        case 'k':
          opt->shm_key = atoi(optarg); break;
        case 'v':
          tmap_progress_set_verbosity(1); break;
        case 'h':
          free(getopt_format);
          return 0;
          break;
        case 'u':
          opt->min_seq_len = atoi(optarg); break;
        case 'U':
          opt->max_seq_len = atoi(optarg); break;
        default:
          // algorithm-specific options
          switch(opt->algo_id) {
            case TMAP_MAP_ALGO_MAP1:
              switch(c) {
                case 'l':
                  opt->seed_length = atoi(optarg); break;
                case 's':
                  opt->seed_max_diff = atoi(optarg); break;
                case 'L':
                  opt->seed2_length = atoi(optarg); break;
                case 'p':
                  if(NULL != strstr(optarg, ".")) opt->max_diff = -1, opt->max_diff_fnr = atof(optarg);
                  else opt->max_diff = atoi(optarg), opt->max_diff_fnr = -1.0;
                  break;
                case 'P':
                  opt->max_err_rate = atof(optarg); break;
                case 'm':
                  if(NULL != strstr(optarg, ".")) opt->max_mm = -1, opt->max_mm_frac = atof(optarg);
                  else opt->max_mm = atoi(optarg), opt->max_mm_frac = -1.0;
                  break;
                case 'o':
                  if(NULL != strstr(optarg, ".")) opt->max_gapo = -1, opt->max_gapo_frac = atof(optarg);
                  else opt->max_gapo = atoi(optarg), opt->max_gapo_frac = -1.0;
                  break;
                case 'e':
                  if(NULL != strstr(optarg, ".")) opt->max_gape = -1, opt->max_gape_frac = atof(optarg);
                  else opt->max_gape = atoi(optarg), opt->max_gape_frac = -1.0;
                  break;
                case 'd':
                  opt->max_cals_del = atoi(optarg); break;
                case 'i':
                  opt->indel_ends_bound = atoi(optarg); break;
                case 'b':
                  opt->max_best_cals = atoi(optarg); break;
                case 'Q':
                  opt->max_entries = atoi(optarg); break;
                default:
                  free(getopt_format);
                  return 0;
              }
              break;
            case TMAP_MAP_ALGO_MAP2:
              switch(c) {
                  /*
                     case 'y': 
                     opt->yita = atof(optarg); break;
                     */
                  /*
                     case 'm': 
                     opt->mask_level = atof(optarg); break;
                     */
                case 'c':
                  opt->length_coef = atof(optarg); break;
                case 'S':
                  opt->max_seed_intv = atoi(optarg); break;
                case 'b':
                  opt->z_best= atoi(optarg); break;
                case 'N':
                  opt->seeds_rev = atoi(optarg); break;
                default: 
                  free(getopt_format);
                  return 0;
              }
              break;
            case TMAP_MAP_ALGO_MAP3:
              switch(c) {
                case 'l':
                  opt->seed_length = atoi(optarg); opt->seed_length_set = 1; break;
                case 'S':
                  opt->max_seed_hits = atoi(optarg); break;
                case 'H':
                  opt->hp_diff = atoi(optarg); break;
                  break;
                default:
                  free(getopt_format);
                  return 0;
              }
              break;
            case TMAP_MAP_ALGO_MAPALL:
              switch(c) {
                case 'I':
                  opt->aln_output_mode_ind = 1; break;
                default:
                  free(getopt_format);
                  return 0;
              }
              break;
            default:
              break;
          }
          break;
      }
  }
  free(getopt_format);
  return 1;
}

static int32_t
tmap_map_opt_file_check_with_null(char *fn1, char *fn2)
{
  if(NULL == fn1 && NULL == fn2) {
      return 0;
  }
  else if((NULL == fn1 && NULL != fn2)
          || (NULL != fn1 && NULL == fn2)) {
      return 1;
  }
  else if(0 != strcmp(fn1, fn2)) {
      return 1;
  }
  return 0;
}

// check that the mapall global options match the algorithm specific global options
#define __tmap_map_opt_check_common1(opt_map_all, opt_map_other) do { \
    if(0 != tmap_map_opt_file_check_with_null((opt_map_other)->fn_fasta, (opt_map_all)->fn_fasta)) { \
        tmap_error("option -f was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if(0 != tmap_map_opt_file_check_with_null((opt_map_other)->fn_reads, (opt_map_all)->fn_reads)) { \
        tmap_error("option -r was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->reads_format != (opt_map_all)->reads_format) { \
        tmap_error("option -F was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->score_match != (opt_map_all)->score_match) { \
        tmap_error("option -A was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->pen_mm != (opt_map_all)->pen_mm) { \
        tmap_error("option -M was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->pen_gapo != (opt_map_all)->pen_gapo) { \
        tmap_error("option -O was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->pen_gape != (opt_map_all)->pen_gape) { \
        tmap_error("option -E was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->fscore != (opt_map_all)->fscore) { \
        tmap_error("option -X was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if(0 != tmap_map_opt_file_check_with_null((opt_map_other)->flow_order, (opt_map_all)->flow_order)) { \
        tmap_error("option -x was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->flow_order_use_sff != (opt_map_all)->flow_order_use_sff) { \
        tmap_error("option -x was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->bw != (opt_map_all)->bw) { \
        tmap_error("option -w was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->softclip_type != (opt_map_all)->softclip_type) { \
        tmap_error("option -g was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->dup_window != (opt_map_all)->dup_window) { \
        tmap_error("option -W was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->max_seed_band != (opt_map_all)->max_seed_band) { \
        tmap_error("option -B was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->score_thr != (opt_map_all)->score_thr) { \
        tmap_error("option -T was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->reads_queue_size != (opt_map_all)->reads_queue_size) { \
        tmap_error("option -q was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->num_threads != (opt_map_all)->num_threads) { \
        tmap_error("option -n was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if(0 != tmap_map_opt_file_check_with_null((opt_map_other)->sam_rg, (opt_map_all)->sam_rg)) { \
        tmap_error("option -R was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->sam_sff_tags != (opt_map_all)->sam_sff_tags) { \
        tmap_error("option -Y was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->input_compr != (opt_map_all)->input_compr) { \
        tmap_error("option -j or -z was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->output_compr != (opt_map_all)->output_compr) { \
        tmap_error("option -J or -Z was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->shm_key != (opt_map_all)->shm_key) { \
        tmap_error("option -s was specified outside of the common options", Exit, CommandLineArgument); \
    } \
} while(0)

void
tmap_map_opt_check(tmap_map_opt_t *opt)
{
  int32_t i;
  // global options
  if(NULL == opt->fn_fasta && 0 == opt->shm_key) {
      tmap_error("option -f or option -s must be specified", Exit, CommandLineArgument);
  }
  else if(NULL != opt->fn_fasta && 0 < opt->shm_key) {
      tmap_error("option -f and option -s may not be specified together", Exit, CommandLineArgument);
  }
  if(NULL == opt->fn_reads && TMAP_READS_FORMAT_UNKNOWN == opt->reads_format) {
      tmap_error("option -F or option -r must be specified", Exit, CommandLineArgument);
  }
  if(TMAP_READS_FORMAT_UNKNOWN == opt->reads_format) {
      tmap_error("the reads format (-r) was unrecognized", Exit, CommandLineArgument);
  }
  tmap_error_cmd_check_int(opt->score_match, 0, INT32_MAX, "-A");
  tmap_error_cmd_check_int(opt->pen_mm, 0, INT32_MAX, "-M");
  tmap_error_cmd_check_int(opt->pen_gapo, 0, INT32_MAX, "-O");
  tmap_error_cmd_check_int(opt->pen_gape, 0, INT32_MAX, "-E");
  tmap_error_cmd_check_int(opt->fscore, 0, INT32_MAX, "-X");
  if(NULL != opt->flow_order) {
      if(0 == strcmp("sff", opt->flow_order) || 0 == strcmp("SFF", opt->flow_order)) {
          if(TMAP_READS_FORMAT_SFF != opt->reads_format) {
              tmap_error("an SFF was not specified (-r) but you want to use the sff flow order (-x)", Exit, CommandLineArgument);
          }
      }
      else {
          tmap_validate_flow_order(opt->flow_order);
      }
  }
  tmap_error_cmd_check_int(opt->bw, 0, INT32_MAX, "-w");
  tmap_error_cmd_check_int(opt->softclip_type, 0, 3, "-g");
  tmap_error_cmd_check_int(opt->dup_window, -1, INT32_MAX, "-W");
  tmap_error_cmd_check_int(opt->max_seed_band, 1, INT32_MAX, "-B");
  tmap_error_cmd_check_int(opt->score_thr, 0, INT32_MAX, "-T");
  if(-1 != opt->reads_queue_size) tmap_error_cmd_check_int(opt->reads_queue_size, 1, INT32_MAX, "-q");
  tmap_error_cmd_check_int(opt->num_threads, 1, INT32_MAX, "-n");
  tmap_error_cmd_check_int(opt->aln_output_mode, 0, 3, "-a");
  if(TMAP_FILE_BZ2_COMPRESSION == opt->output_compr
     && -1 == opt->reads_queue_size) {
      tmap_error("cannot buffer reads with bzip2 output (options \"-q 1 -J\")", Exit, OutOfRange);
  }
  if(-1 != opt->min_seq_len) tmap_error_cmd_check_int(opt->min_seq_len, 1, INT32_MAX, "-u");
  if(-1 != opt->max_seq_len) tmap_error_cmd_check_int(opt->max_seq_len, 1, INT32_MAX, "-U");
  if(-1 != opt->min_seq_len && -1 != opt->max_seq_len &&
     opt->max_seq_len < opt->min_seq_len) {
      tmap_error("The minimum sequence length must be less than the maximum sequence length (-u and -U)", Exit, OutOfRange);
  }

  switch(opt->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      // map1 options
      tmap_error_cmd_check_int((opt->max_diff_fnr < 0) ? opt->max_diff: (int32_t)opt->max_diff_fnr, 0, INT32_MAX, "-p");
      tmap_error_cmd_check_int((int32_t)opt->max_err_rate, 0, INT32_MAX, "-P");
      // this will take care of the case where they are both < 0
      tmap_error_cmd_check_int((opt->max_mm_frac < 0) ? opt->max_mm : (int32_t)opt->max_mm_frac, 0, INT32_MAX, "-m");
      // this will take care of the case where they are both < 0
      tmap_error_cmd_check_int((opt->max_gapo_frac < 0) ? opt->max_gapo : (int32_t)opt->max_gapo_frac, 0, INT32_MAX, "-m");
      // this will take care of the case where they are both < 0
      tmap_error_cmd_check_int((opt->max_gape_frac < 0) ? opt->max_gape : (int32_t)opt->max_gape_frac, 0, INT32_MAX, "-m");
      tmap_error_cmd_check_int(opt->max_cals_del, 1, INT32_MAX, "-d");
      tmap_error_cmd_check_int(opt->indel_ends_bound, 0, INT32_MAX, "-i");
      tmap_error_cmd_check_int(opt->max_best_cals, 0, INT32_MAX, "-b");
      tmap_error_cmd_check_int(opt->max_entries, 1, INT32_MAX, "-Q");
      if(-1 != opt->seed_length) tmap_error_cmd_check_int(opt->seed_length, 1, INT32_MAX, "-l");
      if(-1 != opt->seed2_length) tmap_error_cmd_check_int(opt->seed2_length, 1, INT32_MAX, "-l");
      if(-1 != opt->seed_length && -1 != opt->seed2_length) {
          tmap_error_cmd_check_int(opt->seed_length, 1, opt->seed2_length, "The secondary seed length (-L) must be less than the primary seed length (-l)");
      }
      break;
    case TMAP_MAP_ALGO_MAP2:
      //tmap_error_cmd_check_int(opt->yita, 0, 1, "-y");
      //tmap_error_cmd_check_int(opt->mask_level, 0, 1, "-m");
      tmap_error_cmd_check_int(opt->length_coef, 0, INT32_MAX, "-c");
      tmap_error_cmd_check_int(opt->max_seed_intv, 0, INT32_MAX, "-S");
      tmap_error_cmd_check_int(opt->z_best, 1, INT32_MAX, "-Z");
      tmap_error_cmd_check_int(opt->seeds_rev, 0, INT32_MAX, "-N");
      break;
    case TMAP_MAP_ALGO_MAP3:
      if(-1 != opt->seed_length) tmap_error_cmd_check_int(opt->seed_length, 1, INT32_MAX, "-l");
      tmap_error_cmd_check_int(opt->max_seed_hits, 1, INT32_MAX, "-S");
      tmap_error_cmd_check_int(opt->hp_diff, 0, INT32_MAX, "-H");
      if(0 < opt->hp_diff && TMAP_SEQ_TYPE_SFF != opt->reads_format) tmap_error("-H option must be used with SFF only", Exit, OutOfRange); 
      break;
    case TMAP_MAP_ALGO_MAPALL:
      tmap_error_cmd_check_int(opt->aln_output_mode_ind, 0, 1, "-I");
      if(0 == opt->algos[0] || 0 == opt->num_stages) {
          tmap_error("no algorithms given for stage 1", Exit, CommandLineArgument);
      }
      for(i=0;i<2;i++) {
          // check mapping algorithm specific options
          tmap_map_opt_check(opt->opt_map1[i]);
          tmap_map_opt_check(opt->opt_map2[i]);
          tmap_map_opt_check(opt->opt_map3[i]);

          // check that common values match other opt values
          __tmap_map_opt_check_common1(opt, opt->opt_map1[i]);
          __tmap_map_opt_check_common1(opt, opt->opt_map2[i]);
          __tmap_map_opt_check_common1(opt, opt->opt_map3[i]);
      }
      break;
    default:
      break;
  }
}

void
tmap_map_opt_print(tmap_map_opt_t *opt)
{
  fprintf(stderr, "algo_id=%d\n", opt->algo_id);
  fprintf(stderr, "fn_fasta=%s\n", opt->fn_fasta);
  fprintf(stderr, "fn_reads=%s\n", opt->fn_reads);
  fprintf(stderr, "reads_format=%d\n", opt->reads_format);
  fprintf(stderr, "score_match=%d\n", opt->score_match);
  fprintf(stderr, "pen_mm=%d\n", opt->pen_mm);
  fprintf(stderr, "pen_gapo=%d\n", opt->pen_gapo);
  fprintf(stderr, "pen_gape=%d\n", opt->pen_gape);
  fprintf(stderr, "fscore=%d\n", opt->fscore);
  fprintf(stderr, "flow_order=%s\n", opt->flow_order);
  fprintf(stderr, "flow_order_use_sff=%d\n", opt->flow_order_use_sff);
  fprintf(stderr, "bw=%d\n", opt->bw);
  fprintf(stderr, "softclip_type=%d\n", opt->softclip_type);
  fprintf(stderr, "dup_window=%d\n", opt->dup_window);
  fprintf(stderr, "max_seed_band=%d\n", opt->max_seed_band);
  fprintf(stderr, "score_thr=%d\n", opt->score_thr);
  fprintf(stderr, "reads_queue_size=%d\n", opt->reads_queue_size);
  fprintf(stderr, "num_threads=%d\n", opt->num_threads);
  fprintf(stderr, "aln_output_mode=%d\n", opt->aln_output_mode);
  fprintf(stderr, "sam_rg=%s\n", opt->sam_rg);
  fprintf(stderr, "sam_sff_tags=%d\n", opt->sam_sff_tags);
  fprintf(stderr, "input_compr=%d\n", opt->input_compr);
  fprintf(stderr, "output_compr=%d\n", opt->output_compr);
  fprintf(stderr, "shm_key=%d\n", (int)opt->shm_key);
  fprintf(stderr, "min_seq_len=%d\n", opt->min_seq_len);
  fprintf(stderr, "max_seq_len=%d\n", opt->max_seq_len);
  fprintf(stderr, "seed_length=%d\n", opt->seed_length);
  fprintf(stderr, "seed_length_set=%d\n", opt->seed_length_set);
  fprintf(stderr, "seed_max_diff=%d\n", opt->seed_max_diff);
  fprintf(stderr, "seed2_length=%d\n", opt->seed2_length);
  fprintf(stderr, "max_diff=%d\n", opt->max_diff);
  fprintf(stderr, "max_diff_fnr=%lf\n", opt->max_diff_fnr);
  fprintf(stderr, "max_err_rate=%lf\n", opt->max_err_rate);
  fprintf(stderr, "max_mm=%d\n", opt->max_mm);
  fprintf(stderr, "max_mm_frac=%lf\n", opt->max_mm_frac);
  fprintf(stderr, "max_gapo=%d\n", opt->max_gapo);
  fprintf(stderr, "max_gapo_frac=%lf\n", opt->max_gapo_frac);
  fprintf(stderr, "max_gape=%d\n", opt->max_gape);
  fprintf(stderr, "max_gape_frac=%lf\n", opt->max_gape_frac);
  fprintf(stderr, "max_cals_del=%d\n", opt->max_cals_del);
  fprintf(stderr, "indel_ends_bound=%d\n", opt->indel_ends_bound);
  fprintf(stderr, "max_best_cals=%d\n", opt->max_best_cals);
  fprintf(stderr, "max_entries=%d\n", opt->max_entries);
  fprintf(stderr, "yita=%lf\n", opt->yita);
  fprintf(stderr, "length_coef=%lf\n", opt->length_coef);
  fprintf(stderr, "max_seed_intv=%d\n", opt->max_seed_intv);
  fprintf(stderr, "z_best=%d\n", opt->z_best);
  fprintf(stderr, "seeds_rev=%d\n", opt->seeds_rev);
  fprintf(stderr, "max_seed_hits=%d\n", opt->max_seed_hits);
  fprintf(stderr, "hp_diff=%d\n", opt->hp_diff);
  fprintf(stderr, "aln_output_mode_ind=%d\n", opt->aln_output_mode_ind);
}

void
tmap_map_sam_malloc_aux(tmap_map_sam_t *s, int32_t algo_id)
{
  switch(s->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      s->aux.map1_aux = tmap_calloc(1, sizeof(tmap_map_map1_aux_t), "s->aux.map1_aux");
      break;
    case TMAP_MAP_ALGO_MAP2:
      s->aux.map2_aux = tmap_calloc(1, sizeof(tmap_map_map2_aux_t), "s->aux.map2_aux");
      break;
    case TMAP_MAP_ALGO_MAP3:
      s->aux.map3_aux = tmap_calloc(1, sizeof(tmap_map_map3_aux_t), "s->aux.map3_aux");
      break;
    default:
      break;
  }
}

inline void
tmap_map_sam_destroy_aux(tmap_map_sam_t *s)
{
  switch(s->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      free(s->aux.map1_aux);
      s->aux.map1_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP2:
      free(s->aux.map2_aux);
      s->aux.map2_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP3:
      free(s->aux.map3_aux);
      s->aux.map3_aux = NULL;
      break;
    default:
      break;
  }
}

void
tmap_map_sam_destroy(tmap_map_sam_t *s)
{
  tmap_map_sam_destroy_aux(s);
  free(s->cigar);
  s->cigar = NULL;
  s->n_cigar = 0;
}

tmap_map_sams_t *
tmap_map_sams_init()
{
  tmap_map_sams_t *sams = tmap_calloc(1, sizeof(tmap_map_sams_t), "sams");
  sams->sams = NULL;
  sams->n = 0;
  return sams;
}

void
tmap_map_sams_realloc(tmap_map_sams_t *s, int32_t n)
{
  int32_t i;
  if(n == s->n) return; 
  for(i=n;i<s->n;i++) {
      tmap_map_sam_destroy(&s->sams[i]);
  }
  s->sams = tmap_realloc(s->sams, sizeof(tmap_map_sam_t) * n, "s->sams");
  for(i=s->n;i<n;i++) {
      // nullify
      s->sams[i].algo_id = TMAP_MAP_ALGO_NONE;
      s->sams[i].n_cigar = 0;
      s->sams[i].cigar = NULL;
      s->sams[i].aux.map1_aux = NULL;
      s->sams[i].aux.map2_aux = NULL;
      s->sams[i].aux.map3_aux = NULL;
      s->sams[i].ascore = INT32_MIN;
  }
  s->n = n;
}

void
tmap_map_sams_destroy(tmap_map_sams_t *s)
{
  int32_t i;
  if(NULL == s) return;
  for(i=0;i<s->n;i++) {
      tmap_map_sam_destroy(&s->sams[i]);
  }
  free(s->sams);
  free(s);
}

void
tmap_map_sam_copy_and_nullify(tmap_map_sam_t *dest, tmap_map_sam_t *src)
{
  (*dest) = (*src);
  src->cigar = NULL;
  switch(src->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      src->aux.map1_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP2:
      src->aux.map2_aux = NULL;
      break;
    case TMAP_MAP_ALGO_MAP3:
      src->aux.map3_aux = NULL;
      break;
    default:
      break;
  }
}

static void
tmap_map_sam_print(tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_map_sam_t *sam, int32_t sam_sff_tags)
{
  if(NULL == sam) { // unmapped
      tmap_sam_print_unmapped(tmap_file_stdout, seq, sam_sff_tags);
  }
  else {
      // Note: samtools does not like this value
      if(INT32_MIN == sam->score_subo) {
          sam->score_subo++;
      }
      switch(sam->algo_id) {
        case TMAP_MAP_ALGO_MAP1:
          tmap_sam_print_mapped(tmap_file_stdout, seq, sam_sff_tags, refseq, 
                                sam->strand, sam->seqid, sam->pos,
                                sam->mapq, sam->cigar, sam->n_cigar,
                                sam->score, sam->ascore, sam->algo_id, sam->algo_stage, "");
          break;
        case TMAP_MAP_ALGO_MAP2:
          if(0 < sam->aux.map2_aux->XI) {
              tmap_sam_print_mapped(tmap_file_stdout, seq, sam_sff_tags, refseq, 
                                    sam->strand, sam->seqid, sam->pos,
                                    sam->mapq, sam->cigar, sam->n_cigar,
                                    sam->score, sam->ascore, sam->algo_id, sam->algo_stage, 
                                    "\tXS:i:%d\tXT:i:%d\t\tXF:i:%d\tXE:i:%d\tXI:i:%d",
                                    sam->score_subo,
                                    sam->n_seeds,
                                    sam->aux.map2_aux->XF, sam->aux.map2_aux->XE, 
                                    sam->aux.map2_aux->XI);
          }
          else {
              tmap_sam_print_mapped(tmap_file_stdout, seq, sam_sff_tags, refseq, 
                                    sam->strand, sam->seqid, sam->pos,
                                    sam->mapq, sam->cigar, sam->n_cigar,
                                    sam->score, sam->ascore, sam->algo_id, sam->algo_stage, 
                                    "\tXS:i:%d\tXT:i:%d\tXF:i:%d\tXE:i:%d",
                                    sam->score_subo,
                                    sam->n_seeds,
                                    sam->aux.map2_aux->XF, sam->aux.map2_aux->XE);
          }
          break;
        case TMAP_MAP_ALGO_MAP3:
          tmap_sam_print_mapped(tmap_file_stdout, seq, sam_sff_tags, refseq, 
                                sam->strand, sam->seqid, sam->pos,
                                sam->mapq, sam->cigar, sam->n_cigar,
                                sam->score, sam->ascore, sam->algo_id, sam->algo_stage, 
                                "\tXS:i:%d\tXT:i:%d",
                                sam->score_subo,
                                sam->n_seeds);
          break;
      }
  }
}

void 
tmap_map_sams_print(tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_map_sams_t *sams, int32_t sam_sff_tags) 
{
  int32_t i;
  if(0 < sams->n) {
      if(1 < sams->n) { // sort by alignment score
          tmap_sort_introsort(tmap_map_sam_sort_score, sams->n, sams->sams);
      }
      for(i=0;i<sams->n;i++) {
          tmap_map_sam_print(seq, refseq, &sams->sams[i], sam_sff_tags);
      }
  }
  else {
      tmap_map_sam_print(seq, refseq, NULL, sam_sff_tags);
  }
}

void
tmap_map_sams_filter1(tmap_map_sams_t *sams, int32_t aln_output_mode, int32_t algo_id)
{
  int32_t i, j, k;
  int32_t n_best = 0;
  int32_t best_score, cur_score;
  int32_t best_subo;

  if(sams->n <= 1) {
      return;
  }

  for(i=j=0;i<sams->n;i++) {
      if(TMAP_MAP_ALGO_NONE == algo_id
         || sams->sams[i].algo_id == algo_id) {
          j++;
      }
  }
  if(j <= 1) {
      return;
  }

  best_score = best_subo = INT32_MIN;
  n_best = 0;
  for(i=0;i<sams->n;i++) {
      if(TMAP_MAP_ALGO_NONE == algo_id
         || sams->sams[i].algo_id == algo_id) {
          cur_score = sams->sams[i].score;
          if(best_score < cur_score) {
              if(0 < n_best) {
                  best_subo = best_score;
              }
              best_score = cur_score;
              n_best = 1;
          }
          else if(!(cur_score < best_score)) { // equal
              best_subo = best_score; // more than one mapping
              n_best++;
          }
          else if(best_subo < cur_score) {
              best_subo = cur_score;
          }
          // check sub-optimal
          if(TMAP_MAP_ALGO_MAP2 == sams->sams[i].algo_id
             || TMAP_MAP_ALGO_MAP3 == sams->sams[i].algo_id) {
              cur_score = sams->sams[i].score_subo;
              if(best_subo < cur_score) {
                  best_subo = cur_score;
              }
          }

      }
  }

  // adjust mapping qualities
  if(1 < n_best) {
      for(i=0;i<sams->n;i++) {
          if(TMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              sams->sams[i].mapq = 0;
          }
      }
  }
  else {
      for(i=0;i<sams->n;i++) {
          if(TMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              cur_score = sams->sams[i].score;
              if(cur_score < best_score) { // not the best
                  sams->sams[i].mapq = 0;
              }
          }
      }
  }

  // adjust suboptimal
  if(TMAP_MAP_ALGO_NONE == algo_id) {
      for(i=0;i<sams->n;i++) {
          sams->sams[i].score_subo = best_subo;
      }
  }

  if(TMAP_MAP_UTIL_ALN_MODE_ALL == aln_output_mode) {
      return;
  }

  // copy to the front
  if(n_best < sams->n) {
      for(i=j=0;i<sams->n;i++) {
          if(TMAP_MAP_ALGO_NONE == algo_id
             || sams->sams[i].algo_id == algo_id) {
              cur_score = sams->sams[i].score;
              if(cur_score < best_score) { // not the best
                  tmap_map_sam_destroy(&sams->sams[i]);
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                  }
                  j++;
              }
          }
          else {
              if(j < i) { // copy if we are not on the same index
                  tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
              }
              j++;
          }
      }
      // reallocate
      tmap_map_sams_realloc(sams, j);
  }

  if(TMAP_MAP_UTIL_ALN_MODE_UNIQ_BEST == aln_output_mode) {
      if(1 < n_best) { // there can only be one
          if(TMAP_MAP_ALGO_NONE == algo_id) {
              tmap_map_sams_realloc(sams, 0);
          }
          else {
              // get rid of all of them
              for(i=j=0;i<sams->n;i++) {
                  if(sams->sams[i].algo_id == algo_id) {
                      tmap_map_sam_destroy(&sams->sams[i]);
                  }
                  else {
                      if(j < i) { // copy if we are not on the same index
                          tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                      }
                      j++;
                  }
              }
              tmap_map_sams_realloc(sams, j);
          }
      }
  }
  else if(TMAP_MAP_UTIL_ALN_MODE_RAND_BEST == aln_output_mode) { // get a random
      int32_t r = (int32_t)(drand48() * n_best);

      // keep the rth one
      if(TMAP_MAP_ALGO_NONE == algo_id) {
          if(0 != r) {
              tmap_map_sam_destroy(&sams->sams[0]);
              tmap_map_sam_copy_and_nullify(&sams->sams[0], &sams->sams[r]);
          }
          // reallocate
          tmap_map_sams_realloc(sams, 1);
      }
      else {
          // keep the rth one
          for(i=j=k=0;i<sams->n;i++) {
              if(sams->sams[i].algo_id == algo_id) {
                  if(k == r) { // keep
                      if(j < i) { // copy if we are not on the same index
                          tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                      }
                      j++;
                  }
                  else { // free
                      tmap_map_sam_destroy(&sams->sams[i]);
                  }
                  k++;
              }
              else {
                  if(j < i) { // copy if we are not on the same index
                      tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[i]);
                  }
                  j++;
              }
          }
          tmap_map_sams_realloc(sams, j);
      }
  }
  else if(TMAP_MAP_UTIL_ALN_MODE_ALL_BEST == aln_output_mode) {
      // do nothing
  }
  else {
      tmap_error("bug encountered", Exit, OutOfRange);
  }
}

void
tmap_map_sams_filter(tmap_map_sams_t *sams, int32_t aln_output_mode)
{
  tmap_map_sams_filter1(sams, aln_output_mode, TMAP_MAP_ALGO_NONE);
}

void
tmap_map_util_remove_duplicates(tmap_map_sams_t *sams, int32_t dup_window)
{
  int32_t i, next_i, j, k, end, best_score_i, best_score_n, best_score_subo;

  if(dup_window < 0 || sams->n <= 0) {
      return;
  }

  // sort
  tmap_sort_introsort(tmap_map_sam_sort_coord, sams->n, sams->sams);

  // remove duplicates within a window
  for(i=j=0;i<sams->n;) {

      // get the change
      end = best_score_i = i;
      best_score_n = 0;
      best_score_subo = INT32_MIN;
      while(end+1 < sams->n) {
          if(sams->sams[end].seqid == sams->sams[end+1].seqid
             && sams->sams[end].strand == sams->sams[end+1].strand
             && fabs(sams->sams[end].pos - sams->sams[end+1].pos) <= dup_window) {
              // track the best scoring
              if(sams->sams[best_score_i].score == sams->sams[end+1].score) {
                  best_score_i = end+1;
                  best_score_n++;
              }
              else if(sams->sams[best_score_i].score < sams->sams[end+1].score) {
                  best_score_i = end+1;
                  best_score_n = 1;
              }
              if(best_score_subo < sams->sams[end+1].score_subo) {
                  best_score_subo = sams->sams[end+1].score_subo;
              }
              end++;
          }
          else {
              break;
          }
      }
      next_i = end+1;

      // randomize the best scoring      
      if(1 < best_score_n) {
          k = (int32_t)(best_score_n * drand48()); // make this zero-based 
          best_score_n = 0; // make this one-based
          end = i;
          while(best_score_n <= k) { // this assumes we know there are at least "best_score
              if(sams->sams[best_score_i].score == sams->sams[end].score) {
                  best_score_i = end;
                  best_score_n++;
              }
              end++;
          }
      }

      // copy over the best
      if(j != best_score_i) {
          // destroy
          tmap_map_sam_destroy(&sams->sams[j]);
          // nullify
          tmap_map_sam_copy_and_nullify(&sams->sams[j], &sams->sams[best_score_i]);
      }

      // copy over sub-optimal score
      sams->sams[j].score_subo = best_score_subo;

      // next
      i = next_i;
      j++;
  }

  // resize
  tmap_map_sams_realloc(sams, j);
}

inline void
tmap_map_util_mapq(tmap_map_sams_t *sams, int32_t seq_len, tmap_map_opt_t *opt)
{
  int32_t i, best_i;
  int32_t n_best = 0, n_best_subo = 0;
  int32_t best_score, cur_score, best_subo, best_subo2;
  int32_t mapq;
  int32_t stage = -1;
  int32_t algo_id = TMAP_MAP_ALGO_NONE; 

  // estimate mapping quality TODO: this needs to be refined
  best_i = 0;
  best_score = INT32_MIN;
  best_subo = best_subo2 = opt->score_thr;
  n_best = n_best_subo = 0;
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      if(best_score < cur_score) {
          // save sub-optimal
          best_subo = best_score;
          n_best_subo = n_best;
          // update
          best_score = cur_score;
          n_best = 1;
          best_i = i;
          stage = (algo_id == TMAP_MAP_ALGO_NONE) ? sams->sams[i].algo_stage-1 : -1;
          algo_id = (algo_id == TMAP_MAP_ALGO_NONE) ? sams->sams[i].algo_id : -1;
      }
      else if(cur_score == best_score) { // qual
          n_best++;
      }
      else {
          if(best_subo < cur_score) {
              best_subo = cur_score;
              n_best_subo = 1;
          }
          else if(best_subo == cur_score) {
              n_best_subo++;
          }
      }
      // get the best subo-optimal score
      cur_score = sams->sams[i].score_subo;
      if(INT32_MIN == cur_score) {
          // ignore
      }
      else if(best_subo < cur_score) {
          best_subo2 = cur_score;
      }
  }
  if(best_subo < best_subo2) best_subo = best_subo2;
  if(1 < n_best || best_score <= best_subo) {
      mapq = 0;
  }
  else {
      if(0 == n_best_subo) {
          n_best_subo = 1;
          best_subo = opt->score_thr; 
      }
      /*
         fprintf(stderr, "n_best=%d n_best_subo=%d\n",
         n_best, n_best_subo);
         fprintf(stderr, "best_score=%d best_subo=%d\n",
         best_score, best_subo);
         */
      // Note: this is the old calculationg, based on BWA-long
      //mapq = (int32_t)((n_best / (1.0 * n_best_subo)) * (best_score - best_subo) * (250.0 / best_score + 0.03 / opt->score_match) + .499);
      //
      double sf = 0.4; // initial scaling factor.  Note: 250 * sf is the maximum mapping quality.
      sf *= 250.0 / ((double)opt->score_match * seq_len); // scale based on the best possible alignment score
      sf *= (n_best / (1.0 * n_best_subo)); // scale based on number of sub-optimal mappings
      sf *= (double)(best_score - best_subo + 1 ); // scale based on distance to the sub-optimal mapping
      //sf *= (seq_len < 10) ? 1.0 : log10(seq_len); // scale based on longer reads having more information content
      mapq = (int32_t)(sf + 0.99999);
      if(mapq > 250) mapq = 250;
      if(mapq <= 0) mapq = 1;
  }
  for(i=0;i<sams->n;i++) {
      cur_score = sams->sams[i].score;
      if(cur_score == best_score) {
          sams->sams[i].mapq = mapq;
      }
      else {
          sams->sams[i].mapq = 0;
      }
  }
}

#define __map_util_gen_ap(par, opt) do { \
    int32_t i; \
    for(i=0;i<25;i++) { \
        (par).matrix[i] = -(opt)->pen_mm; \
    } \
    for(i=0;i<4;i++) { \
        (par).matrix[i*5+i] = (opt)->score_match; \
    } \
    (par).gap_open = (opt)->pen_gapo; (par).gap_ext = (opt)->pen_gape; \
    (par).gap_end = (opt)->pen_gape; \
    (par).row = 5; \
} while(0)


tmap_map_sams_t *
tmap_map_util_sw(tmap_refseq_t *refseq,
                 tmap_map_sams_t *sams, 
                 tmap_seq_t *seq,
                 tmap_map_opt_t *opt)
{
  int32_t i, matrix[25];
  int32_t start, end;
  tmap_map_sams_t *sams_tmp = NULL;
  tmap_sw_param_t par;
  tmap_sw_path_t *path = NULL;
  int32_t path_len, path_mem=0;
  tmap_seq_t *seqs[2] = {NULL, NULL};
  int32_t seq_len=0, target_len, target_mem=0;
  uint8_t *target=NULL;
  int32_t best_subo;

  if(0 == sams->n) {
      return sams;
  }

  // the final mappings will go here 
  sams_tmp = tmap_map_sams_init();
  tmap_map_sams_realloc(sams_tmp, sams->n);

  // scoring matrix
  par.matrix = matrix;
  __map_util_gen_ap(par, opt); // TODO

  // sort by strand/chr/pos/score
  tmap_sort_introsort(tmap_map_sam_sort_coord, sams->n, sams->sams);

  i = start = end = 0;
  best_subo = INT32_MIN;
  while(end < sams->n) {
      uint8_t strand, *query=NULL;
      uint32_t start_pos, end_pos;
      tmap_map_sam_t tmp_sam;

      // get the strand/start/end positions
      strand = sams->sams[end].strand;
      start_pos = sams->sams[end].pos + 1; // one-based
      end_pos = sams->sams[end].pos + sams->sams[end].target_len; // one-based

      // update the query sequence
      if(NULL == seqs[strand]) {
          seq_len = tmap_seq_get_bases(seq)->l;
          seqs[strand] = tmap_seq_clone(seq); // clone
          if(1 == strand) { // reverse compliment
              tmap_seq_reverse_compliment(seqs[strand]);
          }
          tmap_seq_to_int(seqs[strand]);
      }
      query = (uint8_t*)tmap_seq_get_bases(seqs[strand])->s;

      // check if the hits can be banded
      if(end + 1 < sams->n) {
          if(sams->sams[end].strand == sams->sams[end+1].strand 
             && sams->sams[end].seqid == sams->sams[end+1].seqid
             && sams->sams[end+1].pos - sams->sams[end].pos <= opt->max_seed_band) {
              end++;
              if(end_pos < sams->sams[end].pos + sams->sams[end].target_len) {
                  end_pos = sams->sams[end].pos + sams->sams[end].target_len; // one-based
              }
              if(best_subo < sams->sams[end].score_subo) {
                  best_subo = sams->sams[end].score_subo;
              }
              continue; // there may be more to add
          }
      }

      // choose a random one within the window
      if(start == end) {
          tmp_sam = sams->sams[start];
      }
      else {
          int32_t r = (int32_t)(drand48() * (end - start + 1));
          r += start;
          tmp_sam = sams->sams[r];
      }

      // update the best sub-optimal score
      tmp_sam.score_subo = best_subo;

      // one-based
      //fprintf(stderr, "1 start_pos=%d\tend_pos=%d\n", start_pos, end_pos);
      if(start_pos < opt->bw) {
          start_pos = 1;
      }
      else {
          start_pos -= opt->bw - 1;
      }
      end_pos += opt->bw - 1;
      if(refseq->annos[sams->sams[end].seqid].len < end_pos) {
          end_pos = refseq->annos[sams->sams[end].seqid].len; // one-based
      }
      //fprintf(stderr, "2 start_pos=%d\tend_pos=%d\n", start_pos, end_pos);

      // get the target sequence
      target_len = end_pos - start_pos + 1;
      if(target_mem < target_len) { // more memory?
          target_mem = target_len;
          tmap_roundup32(target_mem);
          target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }
      if(target_len != tmap_refseq_subseq(refseq, refseq->annos[sams->sams[end].seqid].offset + start_pos, target_len, target)) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }

      // path memory
      if(path_mem <= target_len + seq_len) { // lengthen the path
          path_mem = target_len + seq_len;
          tmap_roundup32(path_mem);
          path = tmap_realloc(path, sizeof(tmap_sw_path_t)*path_mem, "path");
      }

      // get the band width
      par.band_width = end_pos - start_pos + 1;
      par.band_width += 2 * opt->bw; // add bases to the window

      /*
      int j;
      for(j=0;j<target_len;j++) {
          fputc("ACGTN"[target[j]], stderr);
      }
      fputc('\n', stderr);
      for(j=0;j<seq_len;j++) {
          fputc("ACGTN"[query[j]], stderr);
      }
      fputc('\n', stderr);
      */

      // TODO: we could hash previous SWs, to avoid duplicates
      if(0 < tmap_map_util_sw_aux(&tmp_sam,
                              target, target_len,
                              query, seq_len,
                              sams->sams[end].seqid, start_pos-1,
                              &par, path, &path_len,
                              opt->score_thr, opt->softclip_type, strand)) {
          tmap_map_sam_t *s = &sams_tmp->sams[i];
      
          // shallow copy previous data 
          (*s) = tmp_sam; 

          s->n_seeds = (end - start + 1);

          // update aux data
          tmap_map_sam_malloc_aux(s, s->algo_id);
          switch(s->algo_id) {
            case TMAP_MAP_ALGO_MAP1:
              (*s->aux.map1_aux) = (*tmp_sam.aux.map1_aux);
              break;
            case TMAP_MAP_ALGO_MAP2:
              (*s->aux.map2_aux) = (*tmp_sam.aux.map2_aux);
              break;
            case TMAP_MAP_ALGO_MAP3:
              (*s->aux.map3_aux) = (*tmp_sam.aux.map3_aux);
              break;
            default:
              tmap_error("bug encountered", Exit, OutOfRange);
              break;
          }

          i++;
      }

      // update start/end
      end++;
      start = end;
      best_subo = INT32_MIN;
  }
      
  // realloc
  tmap_map_sams_realloc(sams_tmp, i);

  // free memory
  for(i=0;i<2;i++) {
      if(NULL != seqs[i]) {
          tmap_seq_destroy(seqs[i]);
      }
  }
  tmap_map_sams_destroy(sams);
  free(path);
  free(target);

  return sams_tmp;
}

#define _tmap_map_util_sw_aux_path_adjust(_p, _ql, _tl) do { \
    (_p).i = _tl - (_p).i + 1; \
    (_p).j = _ql - (_p).j + 1; \
} while(0) \

int32_t 
tmap_map_util_sw_aux(tmap_map_sam_t *sam,
                 uint8_t *target, int32_t target_length,
                 uint8_t *query, int32_t query_length,
                 uint32_t seqid, uint32_t pos,
                 tmap_sw_param_t *par, tmap_sw_path_t *path, int32_t *path_len,
                 int32_t score_thr, int32_t softclip_type, int32_t strand)
{
  int32_t i, score, score_subo;

  if(1 == strand) { // reverse, but do not compliment, since we want consistent behavior on the nucleotides
      tmap_reverse_int(query, query_length);
      tmap_reverse_int(target, target_length);
  }

  /*
  fprintf(stderr, "strand=%d\n", strand);
  for(i=0;i<query_length;i++)
    fputc("ACGTN"[query[i]], stderr);
  fputc('\n', stderr);
  for(i=0;i<target_length;i++)
    fputc("ACGTN"[target[i]], stderr);
  fputc('\n', stderr);
  */

  switch(softclip_type) {
    case TMAP_MAP_UTIL_SOFT_CLIP_ALL:
      //fprintf(stderr, "TMAP_MAP_UTIL_SOFT_CLIP_ALL\n");
      score = tmap_sw_clipping_core(target, target_length, query, query_length, par, 1, 1, path, path_len, strand);
      break;
    case TMAP_MAP_UTIL_SOFT_CLIP_LEFT:
      score = tmap_sw_clipping_core(target, target_length, query, query_length, par, 1, 0, path, path_len, strand);
      break;
    case TMAP_MAP_UTIL_SOFT_CLIP_RIGHT:
      score = tmap_sw_clipping_core(target, target_length, query, query_length, par, 0, 1, path, path_len, strand);
      break;
    case TMAP_MAP_UTIL_SOFT_CLIP_NONE:
    default:
      score = tmap_sw_clipping_core(target, target_length, query, query_length, par, 0, 0, path, path_len, strand);
      break;
  }
  score_subo = sam->score_subo;

  if(0 < (*path_len) && score_thr < score) {

      if(1 == strand) { // reverse compliment 
          // adjust path
          // TODO

          // reverse the path and adjust its values
          for(i=0;i<(*path_len)>>1;i++) {
              // reverse
              tmap_sw_path_t tmp= path[i];
              path[i] = path[(*path_len)-i-1];
              path[(*path_len)-i-1] = tmp;
              // adjust 
              _tmap_map_util_sw_aux_path_adjust(path[i], query_length, target_length);
              _tmap_map_util_sw_aux_path_adjust(path[(*path_len)-i-1], query_length, target_length);
          }
          if(1 == ((*path_len) & 1)) {
              _tmap_map_util_sw_aux_path_adjust(path[i], query_length, target_length);
          }
      }

      sam->strand = strand;
      sam->seqid = seqid;
      sam->pos = pos + (path[(*path_len)-1].i-1); // zero-based 
      if(path[(*path_len)-1].ctype == TMAP_SW_FROM_I) {
          sam->pos++;
      }
      sam->score = score;
      sam->score_subo = score_subo;
      sam->cigar = tmap_sw_path2cigar(path, (*path_len), &sam->n_cigar);

      if(0 == sam->n_cigar) {
          tmap_error("bug encountered", Exit, OutOfRange);
      }

      // add soft clipping 
      if(1 < path[(*path_len)-1].j) {
          // soft clip the front of the read
          sam->cigar = tmap_realloc(sam->cigar, sizeof(uint32_t)*(1+sam->n_cigar), "sam->cigar");
          for(i=sam->n_cigar-1;0<=i;i--) { // shift up
              sam->cigar[i+1] = sam->cigar[i];
          }
          TMAP_SW_CIGAR_STORE(sam->cigar[0], BAM_CSOFT_CLIP, path[(*path_len)-1].j-1);
          sam->n_cigar++;
      }
      if(path[0].j < query_length) {
          // soft clip the end of the read
          sam->cigar = tmap_realloc(sam->cigar, sizeof(uint32_t)*(1+sam->n_cigar), "sam->cigar");
          TMAP_SW_CIGAR_STORE(sam->cigar[sam->n_cigar], BAM_CSOFT_CLIP, query_length - path[0].j);
          sam->n_cigar++;
      }
  }
  else {
      sam->score = INT32_MIN;
      sam->cigar = NULL;
      sam->n_cigar = 0;
      (*path_len) = 0;
  }
  if(1 == strand) { // reverse back
      tmap_reverse_int(query, query_length);
      tmap_reverse_int(target, target_length);
  }
  return (0 == sam->n_cigar) ? 0 : 1;
}

// TODO: make sure the "longest" read alignment is found
void
tmap_map_util_fsw(tmap_seq_t *seq, 
                  uint8_t *flow_order, int32_t flow_order_len,
                  tmap_map_sams_t *sams, tmap_refseq_t *refseq,
                  int32_t bw, int32_t softclip_type, int32_t score_thr,
                  int32_t score_match, int32_t pen_mm, int32_t pen_gapo, 
                  int32_t pen_gape, int32_t fscore)
{
  int32_t i, j, k, l;
  uint8_t *target = NULL;
  int32_t target_mem = 0, target_len = 0;

  tmap_fsw_flowseq_t *fseq[2] = {NULL, NULL};
  tmap_fsw_path_t *path = NULL;
  int32_t path_mem = 0, path_len = 0;
  tmap_fsw_param_t param;
  int32_t matrix[25];
  tmap_seq_t *seq_rc = NULL;
  uint8_t *flow_order_rc=NULL;
  int32_t softclip_types[2];

  // generate the alignment parameters
  param.matrix = matrix;
  param.band_width = 0;
  param.offset = TMAP_MAP_UTIL_FSW_OFFSET; // this sets the hp difference
  __tmap_fsw_gen_ap1(param, score_match, pen_mm, pen_gapo, pen_gape, fscore);

  softclip_types[0] = softclip_type;
  softclip_types[1] = __tmap_map_util_reverse_soft_clipping(softclip_type);

  // go through each hit
  for(i=0;i<sams->n;i++) {
      tmap_map_sam_t *s = &sams->sams[i];
      uint32_t ref_start, ref_end;

      // get flow sequence if necessary
      if(NULL == fseq[s->strand]) {
          if(0 == s->strand) {
              fseq[s->strand] = tmap_fsw_seq_to_flowseq(seq, flow_order, flow_order_len);
          }
          else {
              flow_order_rc = tmap_malloc(sizeof(uint8_t) * flow_order_len, "flow_order_rc");
              for(j=0;j<flow_order_len;j++) {
                  flow_order_rc[j] = 3 - flow_order[flow_order_len-j-1];
              }
              seq_rc = tmap_seq_clone(seq);
              tmap_seq_reverse_compliment(seq_rc);
              fseq[s->strand] = tmap_fsw_seq_to_flowseq(seq_rc, flow_order_rc, flow_order_len);
          }
          // HERE
          //tmap_fsw_flowseq_print(tmap_file_stderr, fseq[s->strand]);
      }

      // HERE
      /*
      uint32_t *old_cigar, old_n_cigar;
      old_n_cigar = s->n_cigar;
      old_cigar = tmap_malloc(sizeof(uint32_t)*old_n_cigar, "old_cigar");
      for(j=0;j<s->n_cigar;j++) {
          old_cigar[j] = s->cigar[j];
      }
      */

      param.band_width = 0;
      ref_start = ref_end = s->pos + 1;
      for(j=0;j<s->n_cigar;j++) {
          int32_t op, op_len;

          op = TMAP_SW_CIGAR_OP(s->cigar[j]);
          op_len = TMAP_SW_CIGAR_LENGTH(s->cigar[j]);

          switch(op) {
            case BAM_CMATCH:
              ref_end += op_len;
              break;
            case BAM_CDEL:
              if(param.band_width < op_len) param.band_width += op_len;
              ref_end += op_len;
              break;
            case BAM_CINS:
              if(param.band_width < op_len) param.band_width += op_len;
              break;
            case BAM_CSOFT_CLIP:
              if(0 == j) {
                  if(ref_start <= op_len) {
                      ref_start = 1;
                  }
                  else {
                      ref_start = ref_start - op_len + 1;
                  }
              }
              else ref_end += op_len;
              break;
            default:
              // ignore
              break;
          }
      }

      // check bounds
      if(ref_start < 1) ref_start = 1;
      if(refseq->annos[s->seqid].len < ref_end) {
          ref_end = refseq->annos[s->seqid].len;
      }

      // get the target sequence
      target_len = ref_end - ref_start + 1;
      if(target_mem < target_len) {
          target_mem = target_len;
          tmap_roundup32(target_mem);
          target = tmap_realloc(target, sizeof(uint8_t)*target_mem, "target");
      }
      target_len = tmap_refseq_subseq(refseq, ref_start + refseq->annos[s->seqid].offset, target_len, target);

      // add to the band width
      param.band_width += 2 * bw;

      // make sure we have enough memory for the path
      while(path_mem <= target_len + fseq[s->strand]->num_flows) { // lengthen the path
          path_mem = target_len + fseq[s->strand]->num_flows + 1;
          tmap_roundup32(path_mem);
          path = tmap_realloc(path, sizeof(tmap_fsw_path_t)*path_mem, "path");
      }

      /*
      fprintf(stderr, "ref_start=%d ref_end=%d\n", ref_start, ref_end);
      fprintf(stderr, "base_calls:\n");
      for(j=0;j<fseq[s->strand]->num_flows;j++) {
          for(k=0;k<fseq[s->strand]->base_calls[j];k++) {
              fputc("ACGTN"[fseq[s->strand]->flow_order[j % fseq[s->strand]->flow_order_len]], stderr);
          }
      }
      fputc('\n', stderr);
      fprintf(stderr, "target:\n");
      for(j=0;j<target_len;j++) {
          fputc("ACGTN"[target[j]], stderr);
      }
      fputc('\n', stderr);
      for(j=0;j<fseq[s->strand]->flow_order_len;j++) {
          fputc("ACGTN"[fseq[s->strand]->flow_order[j]], stderr);
      }
      fputc('\n', stderr);
      */

      // re-align
      s->ascore = s->score;
      path_len = path_mem;
      //fprintf(stderr, "old score=%d\n", s->score);
      switch(softclip_types[s->strand]) {
        case TMAP_MAP_UTIL_SOFT_CLIP_ALL:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq[s->strand], &param, 
                                            1, 1, path, &path_len);
          break;
        case TMAP_MAP_UTIL_SOFT_CLIP_LEFT:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq[s->strand], &param, 
                                            1, 0, path, &path_len);
          break;
        case TMAP_MAP_UTIL_SOFT_CLIP_RIGHT:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq[s->strand], &param, 
                                            0, 1, path, &path_len);
          break;
        case TMAP_MAP_UTIL_SOFT_CLIP_NONE:
          s->score = tmap_fsw_clipping_core(target, target_len, fseq[s->strand], &param, 
                                            0, 0, path, &path_len);
          break;
        default:
          tmap_error("soft clipping type was not recognized", Exit, OutOfRange);
          break;
      }
      s->score_subo = INT32_MIN;
      //fprintf(stderr, "new score=%d path_len=%d\n", s->score, path_len);

      if(0 < path_len) { // update

          /*
          for(j=0;j<path_len;j++) {
              fprintf(stderr, "j=%d path[j].i=%d path[j].j=%d\n", j, path[j].i, path[j].j);
          }
          */

          s->score = (int32_t)((s->score + 99.99)/100.0); 
          s->pos = (ref_start-1);
          if(0 <= path[path_len-1].j) { // does not begin with an insertion
              s->pos += (path[path_len-1].j);
          }
          if(refseq->len < s->pos) {
              //fprintf(stderr, "ref_start=%d path[path_len-1].j=%d\n", ref_start, path[path_len-1].j);
              tmap_error("bug encountered", Exit, OutOfRange);
          }
          //fprintf(stderr, "path[path_len-1].i=%d path[0].i=%d num_flows=%d\n", path[path_len-1].i, path[0].i, fseq[s->strand]->num_flows);
          free(s->cigar);
          s->cigar = tmap_fsw_path2cigar(path, path_len, &s->n_cigar, 1);

          if(0 < path[path_len-1].i) { // skipped beginning flows
              // get the number of bases to clip
              for(j=k=0;j<path[path_len-1].i;j++) {
                  k += fseq[s->strand]->base_calls[j];
              }
              if(0 < k) { // bases should be soft-clipped
                  s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
                  for(l=s->n_cigar-1;0<=l;l--) {
                      s->cigar[l+1] = s->cigar[l];
                  }
                  TMAP_SW_CIGAR_STORE(s->cigar[0], BAM_CSOFT_CLIP, k);
                  s->n_cigar++;
              }
          }

          if(path[0].i+1 < fseq[s->strand]->num_flows) { // skipped ending flows
              // get the number of bases to clip 
              for(j=path[0].i+1,k=0;j<fseq[s->strand]->num_flows;j++) {
                  k += fseq[s->strand]->base_calls[j];
              }
              if(0 < k) { // bases should be soft-clipped
                  s->cigar = tmap_realloc(s->cigar, sizeof(uint32_t)*(1 + s->n_cigar), "s->cigar");
                  s->cigar[s->n_cigar] = (k << 4) | 4;
                  s->n_cigar++;
              }
          }

          // HERE
              /*
          int32_t differs = 0;
          if(s->n_cigar != old_n_cigar) {
              differs = 1;
          }
          else {
              for(j=0;j<s->n_cigar;j++) {
                  if(s->cigar[j] != old_cigar[j]) {
                      differs = 1;
                      break;
                  }
              }
          }
          if(1 == differs) {
              for(j=0;j<path_len;j++) {
                  fprintf(stderr, "j=%d path[j].i=%d path[j].j=%d\n", j, path[j].i, path[j].j);
              }
              fprintf(stderr, "NEW: [");
              for(j=0;j<s->n_cigar;j++) {
                  fprintf(stderr, "%d%c", s->cigar[j]>>4, "MIDNSHP"[s->cigar[j]&0xf]);
              }
              fprintf(stderr, "]\n");
              fprintf(stderr, "OLD: [");
              for(j=0;j<old_n_cigar;j++) {
                  fprintf(stderr, "%d%c", old_cigar[j]>>4, "MIDNSHP"[old_cigar[j]&0xf]);
              }
              fprintf(stderr, "]\n\n");
          }
              */
      }

      // HERE
      //free(old_cigar);
  }
  // free
  if(NULL != fseq[0]) {
      tmap_fsw_flowseq_destroy(fseq[0]);
  }
  if(NULL != fseq[1]) {
      free(flow_order_rc);
      tmap_fsw_flowseq_destroy(fseq[1]);
      tmap_seq_destroy(seq_rc);
  }
  free(target);
  free(path);
}
