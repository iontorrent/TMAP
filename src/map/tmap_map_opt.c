/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <getopt.h>
#include "../util/tmap_alloc.h"
#include "../util/tmap_error.h"
#include "../io/tmap_file.h"
#include "../seq/tmap_seq.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_definitions.h"
#include "tmap_map_opt.h"


static struct option tmap_map_opt_long_options[] =
{  
    // global options
    {"fn-fasta", required_argument, 0, 'f'},
    {"fn-reads", required_argument, 0, 'r'},
    {"reads-format", required_argument, 0, 'F'},
    {"fn-sam", required_argument, 0, '0'},
    {"score-match", required_argument, 0, 'A'},
    {"pen-mismatch", required_argument, 0, 'M'},
    {"pen-gap-open", required_argument, 0, 'O'},
    {"pen-gap-extension", required_argument, 0, 'E'},
    {"flow-score", required_argument, 0, 'X'},
    {"flow-order", required_argument, 0, 'x'},
    {"key-seq", required_argument, 0, 't'},
    {"band-width", required_argument, 0, 'w'},
    {"softclip-type", required_argument, 0, 'g'},
    {"softclip-key", no_argument, 0, 'y'},
    {"duplicate-window", required_argument, 0, 'W'},
    {"max-seed-band", required_argument, 0, 'B'},
    {"score-thres", required_argument, 0, 'T'},
    {"reads-queue-size", required_argument, 0, 'q'},
    {"num-threads", required_argument, 0, 'n'},
    {"aln-output-mode", required_argument, 0, 'a'},
    {"sam-rg", required_argument, 0, 'R'},
    {"sam-sff-tags", no_argument, 0, 'Y'},
    {"remove-sff-clipping", no_argument, 0, 'G'},
    {"input-gz", no_argument, 0, 'z'},
#ifndef DISABLE_BZ2
    {"input-bz2", no_argument, 0, 'j'},
#endif
    {"output-gz", no_argument, 0, 'Z'},
#ifndef DISABLE_BZ2
    {"output-bz2", no_argument, 0, 'J'},
#endif
    {"shm-key", required_argument, 0, 'k'},
    {"help", no_argument, 0, 'h'},
    {"verbose", no_argument, 0, 'v'},

    // map1/map3 options
    {"seed-length", required_argument, 0, 'l'},

    // map1 options
    {"seed-max-diff", required_argument, 0, 's'},
    {"seed2-length", required_argument, 0, 'L'},
    {"max-diff", required_argument, 0, 'p'},
    {"max-diff-fnr", required_argument, 0, 'p'},
    {"max-err-rate", required_argument, 0, 'P'},
    {"max-mismatch", required_argument, 0, 'm'},
    {"max-mismatch-frac", required_argument, 0, 'm'},
    {"max-gap-open", required_argument, 0, 'o'},
    {"max-gap-open-frac", required_argument, 0, 'o'},
    {"max-gap-extension", required_argument, 0, 'e'},
    {"max-gap-extension-frac", required_argument, 0, 'e'},
    {"max-cals-del", required_argument, 0, 'd'},
    {"indel-ends-bound", required_argument, 0, 'i'},
    {"max-best-cals", required_argument, 0, 'b'},
    {"max-nodes", required_argument, 0, 'Q'},

    // map2 options
    {"length-coef", required_argument, 0, 'c'},
    {"max-seed-intv", required_argument, 0, 'S'},
    {"z-best", required_argument, 0, 'b'},
    {"seeds-rev", required_argument, 0, 'N'},

    // map3 options
    {"max-seed-hits", required_argument, 0, 'S'},
    {"hp-diff", required_argument, 0, 'H'},
    {"hit-frac", required_argument, 0, 'V'},

    // mapvsw options
    // None

    // mapall options
    {"mapall-aln-output-mode-independent", no_argument, 0, 'I'},
    {"mapall-score-thr", required_argument, 0, 'C'},
    {"mapall-mapq-thr", required_argument, 0, 'D'},
    {"mapall-keep-all", no_argument, 0, 'K'},
};

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
  opt->fn_fasta = NULL;
  opt->fn_reads = NULL;
  opt->fn_reads_num = 0;
  opt->reads_format = TMAP_READS_FORMAT_UNKNOWN;
  opt->fn_sam = NULL;
  opt->score_match = TMAP_MAP_OPT_SCORE_MATCH;
  opt->pen_mm = TMAP_MAP_OPT_PEN_MM;
  opt->pen_gapo = TMAP_MAP_OPT_PEN_GAPO;
  opt->pen_gape = TMAP_MAP_OPT_PEN_GAPE;
  opt->fscore = TMAP_MAP_OPT_FSCORE;
  opt->flow_order = NULL;
  opt->flow_order_use_sff = 0;
  opt->key_seq = NULL;
  opt->key_seq_use_sff = 0;
  opt->bw = 50; 
  opt->softclip_type = TMAP_MAP_OPT_SOFT_CLIP_RIGHT;
  opt->softclip_key = 0;
  opt->remove_sff_clipping = 1;
  opt->dup_window = 128;
  opt->max_seed_band = 15;
  opt->score_thr = 8;
  opt->reads_queue_size = 262144;
  opt->num_threads = 1;
  opt->aln_output_mode = TMAP_MAP_OPT_ALN_MODE_RAND_BEST;
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
      opt->seed2_length = 48;
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
      opt->hit_frac = 0.25;
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      // mapvsw
      break;
    case TMAP_MAP_ALGO_MAPALL:
      // mapall
      opt->aln_output_mode_ind = 0;
      opt->mapall_score_thr = 8;
      opt->mapall_mapq_thr = 23; // 0.5% error
      opt->mapall_keep_all = 1;
      for(i=0;i<2;i++) {
          opt->algos[i] = 0;
          opt->opt_map1[i] = tmap_map_opt_init(TMAP_MAP_ALGO_MAP1);
          opt->opt_map2[i] = tmap_map_opt_init(TMAP_MAP_ALGO_MAP2);
          opt->opt_map3[i] = tmap_map_opt_init(TMAP_MAP_ALGO_MAP3);
          opt->opt_map_vsw[i] = tmap_map_opt_init(TMAP_MAP_ALGO_MAPVSW);
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
  for(i=0;i<opt->fn_reads_num;i++) {
      free(opt->fn_reads[i]); 
  }
  free(opt->fn_reads);
  free(opt->fn_sam);
  free(opt->sam_rg);
  free(opt->flow_order);
  free(opt->key_seq);

  switch(opt->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
    case TMAP_MAP_ALGO_MAP2:
    case TMAP_MAP_ALGO_MAP3:
    case TMAP_MAP_ALGO_MAPVSW:
      break;
    case TMAP_MAP_ALGO_MAPALL:
      // mapall
      for(i=0;i<2;i++) {
          opt->algos[i] = 0;
          tmap_map_opt_destroy(opt->opt_map1[i]);
          tmap_map_opt_destroy(opt->opt_map2[i]);
          tmap_map_opt_destroy(opt->opt_map3[i]);
          tmap_map_opt_destroy(opt->opt_map_vsw[i]);
      }
      break;
    default:
      break;
  }

  free(opt);
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
    tmap_file_fprintf(tmap_file_stderr, "         -V INT      the fraction of seed positions that are under the maximum (-S) [%.3lf]\n", (_opt)->hit_frac); \
    tmap_file_fprintf(tmap_file_stderr, "         -u INT      the minimum sequence length to examine [%d]\n", (_opt)->min_seq_len); \
    tmap_file_fprintf(tmap_file_stderr, "         -U INT      the maximum sequence length to examine [%d]\n", (_opt)->max_seq_len); \
    tmap_file_fprintf(tmap_file_stderr, "\n"); \
    } while(0)

#define __tmap_map_opt_usage_map_vsw(_opt, _stage) do { \
    if(_stage < 0) tmap_file_fprintf(tmap_file_stderr, "%s options (optional):\n", tmap_algo_id_to_name(TMAP_MAP_ALGO_MAPVSW)); \
    else tmap_file_fprintf(tmap_file_stderr, "%s stage %d options (optional):\n", tmap_algo_id_to_name(TMAP_MAP_ALGO_MAPVSW), _stage); \
    tmap_file_fprintf(tmap_file_stderr, "         None\n"); \
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
  tmap_file_fprintf(tmap_file_stderr, "         -r FILE     the reads file name [");
  if(0 == opt->fn_reads_num) tmap_file_fprintf(tmap_file_stderr, "stdin");
  else {
      for(i=0;i<opt->fn_reads_num;i++) {
          if(0 < i) tmap_file_fprintf(tmap_file_stderr, ",");
          tmap_file_fprintf(tmap_file_stderr, "%s", opt->fn_reads[i]);
      }
  }
  tmap_file_fprintf(tmap_file_stderr, "]\n");
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "global options (optional):\n");
#ifdef HAVE_SAMTOOLS
  tmap_file_fprintf(tmap_file_stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff|sam|bam) [%s]\n", reads_format);
#else
  tmap_file_fprintf(tmap_file_stderr, "         -F STRING   the reads file format (fastq|fq|fasta|fa|sff) [%s]\n", reads_format);
#endif
  tmap_file_fprintf(tmap_file_stderr, "         -0 FILE     the SAM file name [%s]\n", (NULL == opt->fn_sam) ? "stdout" : opt->fn_sam);
  tmap_file_fprintf(tmap_file_stderr, "         -A INT      score for a match [%d]\n", opt->score_match);
  tmap_file_fprintf(tmap_file_stderr, "         -M INT      the mismatch penalty [%d]\n", opt->pen_mm);
  tmap_file_fprintf(tmap_file_stderr, "         -O INT      the indel start penalty [%d]\n", opt->pen_gapo);
  tmap_file_fprintf(tmap_file_stderr, "         -E INT      the indel extend penalty [%d]\n", opt->pen_gape);
  tmap_file_fprintf(tmap_file_stderr, "         -X INT      the flow score penalty [%d]\n", opt->fscore);
  tmap_file_fprintf(tmap_file_stderr, "         -x STRING   the flow order ([ACGT]{4+} or \"sff\") [%s]\n",
                    (NULL == opt->flow_order) ? "not using" : opt->flow_order);
  tmap_file_fprintf(tmap_file_stderr, "         -t STRING   the key sequence ([ACGT]{4+} or \"sff\") [%s]\n",
                    (NULL == opt->key_seq) ? "not using" : opt->key_seq);
  tmap_file_fprintf(tmap_file_stderr, "         -w INT      the band width [%d]\n", opt->bw);
  tmap_file_fprintf(tmap_file_stderr, "         -g INT      the soft-clipping type [%d]\n", opt->softclip_type);
  tmap_file_fprintf(tmap_file_stderr, "                             0 - allow on the right and left portions of the read\n");
  tmap_file_fprintf(tmap_file_stderr, "                             1 - allow on the left portion of the read\n");
  tmap_file_fprintf(tmap_file_stderr, "                             2 - allow on the right portion of the read\n");
  tmap_file_fprintf(tmap_file_stderr, "                             3 - do not allow soft-clipping\n");
  tmap_file_fprintf(tmap_file_stderr, "         -y          soft clip only the last base of the key [%s]\n",
                    (1 == opt->softclip_key) ? "true" : "false");
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
  tmap_file_fprintf(tmap_file_stderr, "         -G          do not remove SFF clipping [%d]\n", opt->remove_sff_clipping);
  tmap_file_fprintf(tmap_file_stderr, "         -Y          include SFF specific SAM tags [%s]\n",
                    (1 == opt->sam_sff_tags) ? "true" : "false");

#ifndef DISABLE_BZ2
  tmap_file_fprintf(tmap_file_stderr, "         -z/-j       the input is gz/bz2 compressed (gzip/bzip2)");
  __tmap_map_print_compression(opt->input_compr);
  tmap_file_fprintf(tmap_file_stderr, "         -Z/-J       the output is gz/bz2 compressed (gzip/bzip2)");
  __tmap_map_print_compression(opt->output_compr);
#else
  tmap_file_fprintf(tmap_file_stderr, "         -z       the input is gz compressed (gzip)");
  __tmap_map_print_compression(opt->input_compr);
  tmap_file_fprintf(tmap_file_stderr, "         -Z       the output is gz compressed (gzip)");
  __tmap_map_print_compression(opt->output_compr);
#endif
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
    case TMAP_MAP_ALGO_MAPVSW:
      __tmap_map_opt_usage_map_vsw(opt, -1);
      break;
    case TMAP_MAP_ALGO_MAPALL:
      for(i=0;i<2;i++) {
          __tmap_map_opt_usage_map1(opt->opt_map1[i], i+1);
          __tmap_map_opt_usage_map2(opt->opt_map2[i], i+1);
          __tmap_map_opt_usage_map3(opt->opt_map3[i], i+1);
          __tmap_map_opt_usage_map_vsw(opt->opt_map_vsw[i], i+1);
      }
      // mapall
      tmap_file_fprintf(tmap_file_stderr, "%s options (optional):\n", tmap_algo_id_to_name(opt->algo_id));
      tmap_file_fprintf(tmap_file_stderr, "         -I          apply the output filter (-a) and duplicate removal (-W) for each algorithm separately [%s]\n",
                        (1 == opt->aln_output_mode_ind) ? "true" : "false");
      tmap_file_fprintf(tmap_file_stderr, "         -C INT      score threshold for stage one divided by the match score [%d]\n", opt->mapall_score_thr);
      tmap_file_fprintf(tmap_file_stderr, "         -D INT      mapping quality threshold for stage one divided by the match score [%d]\n", opt->mapall_mapq_thr);
      tmap_file_fprintf(tmap_file_stderr, "         -K          do not keep mappings from the first stage for the next stage [%s]\n", 
                        (0 == opt->mapall_keep_all) ? "true" : "false");
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
  int i, c, option_index;
  char *getopt_format = NULL;

  opt->argc = argc; opt->argv = argv;
  switch(opt->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      getopt_format = tmap_strdup("f:r:F:0:A:M:O:E:X:x:t:w:g:yW:B:T:q:n:a:R:YGjzJZk:vhl:s:L:p:P:m:o:e:d:i:b:Q:u:U:");
      break;
    case TMAP_MAP_ALGO_MAP2:
      getopt_format = tmap_strdup("f:r:F:0:A:M:O:E:X:x:t:w:g:yW:B:T:q:n:a:R:YGjzJZk:vhc:S:b:N:u:U:");
      break;
    case TMAP_MAP_ALGO_MAP3:
      getopt_format = tmap_strdup("f:r:F:0:A:M:O:E:X:x:t:w:g:yW:B:T:q:n:a:R:YGjzJZk:vhl:S:H:V:u:U:");
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      getopt_format = tmap_strdup("f:r:F:0:A:M:O:E:X:x:t:w:g:yW:B:T:q:n:a:R:YGjzJZk:vhu:U:");
      break;
    case TMAP_MAP_ALGO_MAPALL:
      getopt_format = tmap_strdup("f:r:F:0:A:M:O:E:X:x:t:w:g:yW:B:T:q:n:a:R:YGjzJZk:vhIC:D:K:u:U:");
      break;
    default:
      tmap_error("unrecognized algorithm", Exit, OutOfRange);
      break;
  }

  // global options
  // Note: possible memory leaks if the same option (besides -R) are specified twice
  while((c = getopt_long(argc, argv, getopt_format, tmap_map_opt_long_options, &option_index)) >= 0) {
      switch(c) { 
        case 'f': 
          opt->fn_fasta = tmap_strdup(optarg); break;
        case 'r':
          opt->fn_reads_num++;
          opt->fn_reads = tmap_realloc(opt->fn_reads, sizeof(char*) * opt->fn_reads_num, "opt->fn_reads");
          opt->fn_reads[opt->fn_reads_num-1] = tmap_strdup(optarg); 
          tmap_get_reads_file_format_from_fn_int(opt->fn_reads[opt->fn_reads_num-1], &opt->reads_format, &opt->input_compr);
          break;
        case 'F':
          opt->reads_format = tmap_get_reads_file_format_int(optarg); break;
        case '0':
          opt->fn_sam = tmap_strdup(optarg); break;
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
        case 't':
          opt->key_seq = tmap_strdup(optarg);
          if(0 == strcmp("sff", opt->key_seq) || 0 == strcmp("SFF", opt->key_seq)) {
              opt->key_seq_use_sff = 1;
          }
          else {
              opt->key_seq_use_sff = 0;
          }
          break;
        case 'w':
          opt->bw = atoi(optarg); break;
        case 'g':
          opt->softclip_type = atoi(optarg); break;
        case 'y':
          opt->softclip_key = 1; break;
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
        case 'G':
          opt->remove_sff_clipping = 0; break;
        case 'j':
          opt->input_compr = TMAP_FILE_BZ2_COMPRESSION;
          for(i=0;i<opt->fn_reads_num;i++) {
              tmap_get_reads_file_format_from_fn_int(opt->fn_reads[i], &opt->reads_format, &opt->input_compr);
          }
          break;
        case 'z':
          opt->input_compr = TMAP_FILE_GZ_COMPRESSION;
          for(i=0;i<opt->fn_reads_num;i++) {
              tmap_get_reads_file_format_from_fn_int(opt->fn_reads[i], &opt->reads_format, &opt->input_compr);
          }
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
                  opt->seed_length = atoi(optarg); 
                  if(0 < opt->seed_length) opt->seed_length_set = 1; 
                  break;
                case 'S':
                  opt->max_seed_hits = atoi(optarg); break;
                case 'H':
                  opt->hp_diff = atoi(optarg); break;
                case 'V':
                  opt->hit_frac = atof(optarg); break;
                default:
                  free(getopt_format);
                  return 0;
              }
              break;
            case TMAP_MAP_ALGO_MAPVSW:
              switch(c) {
                default:
                  free(getopt_format);
                  return 0;
              }
              break;
            case TMAP_MAP_ALGO_MAPALL:
              switch(c) {
                case 'I':
                  opt->aln_output_mode_ind = 1; break;
                case 'C':
                  opt->mapall_score_thr = atoi(optarg); break;
                case 'D':
                  opt->mapall_mapq_thr = atoi(optarg); break;
                case 'K':
                  opt->mapall_keep_all = atoi(optarg); break;
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
    int32_t _i; \
    if(0 != tmap_map_opt_file_check_with_null((opt_map_other)->fn_fasta, (opt_map_all)->fn_fasta)) { \
        tmap_error("option -f was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->fn_reads_num != (opt_map_all)->fn_reads_num) { \
        tmap_error("option -r was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    for(_i=0;_i<(opt_map_other)->fn_reads_num;_i++) { \
        if(0 != tmap_map_opt_file_check_with_null((opt_map_other)->fn_reads[_i], (opt_map_all)->fn_reads[_i])) { \
            tmap_error("option -r was specified outside of the common options", Exit, CommandLineArgument); \
        } \
    } \
    if(0 != tmap_map_opt_file_check_with_null((opt_map_other)->fn_sam, (opt_map_all)->fn_sam)) { \
        tmap_error("option -0 was specified outside of the common options", Exit, CommandLineArgument); \
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
    if(0 != tmap_map_opt_file_check_with_null((opt_map_other)->key_seq, (opt_map_all)->key_seq)) { \
        tmap_error("option -t was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->key_seq_use_sff != (opt_map_all)->key_seq_use_sff) { \
        tmap_error("option -t was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->bw != (opt_map_all)->bw) { \
        tmap_error("option -w was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->softclip_type != (opt_map_all)->softclip_type) { \
        tmap_error("option -g was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->softclip_key != (opt_map_all)->softclip_key) { \
        tmap_error("option -y was specified outside of the common options", Exit, CommandLineArgument); \
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
    if((opt_map_other)->remove_sff_clipping != (opt_map_all)->remove_sff_clipping) { \
        tmap_error("option -G was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->input_compr != (opt_map_all)->input_compr) { \
        tmap_error("option -j or -z was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->output_compr != (opt_map_all)->output_compr) { \
        tmap_error("option -J or -Z was specified outside of the common options", Exit, CommandLineArgument); \
    } \
    if((opt_map_other)->shm_key != (opt_map_all)->shm_key) { \
        tmap_error("option -k was specified outside of the common options", Exit, CommandLineArgument); \
    } \
} while(0)

void
tmap_map_opt_check(tmap_map_opt_t *opt)
{
  int32_t i;
  // global options
  if(NULL == opt->fn_fasta && 0 == opt->shm_key) {
      tmap_error("option -f or option -k must be specified", Exit, CommandLineArgument);
  }
  else if(NULL != opt->fn_fasta && 0 < opt->shm_key) {
      tmap_error("option -f and option -k may not be specified together", Exit, CommandLineArgument);
  }
  if(0 == opt->fn_reads_num && TMAP_READS_FORMAT_UNKNOWN == opt->reads_format) {
      tmap_error("option -F or option -r must be specified", Exit, CommandLineArgument);
  }
  else if(1 < opt->fn_reads_num) {
      if(1 == opt->sam_sff_tags) {
          tmap_error("options -1 and -2 cannot be used with -Y", Exit, CommandLineArgument);
      }
      else if(1 == opt->flow_order_use_sff) {
          tmap_error("options -1 and -2 cannot be used with -x", Exit, CommandLineArgument);
      }
      else if(1 == opt->key_seq_use_sff) {
          tmap_error("options -1 and -2 cannot be used with -t", Exit, CommandLineArgument);
      }
      else if(NULL != opt->flow_order) {
          tmap_error("options -1 and -2 cannot be used with -x", Exit, CommandLineArgument);
      }
      else if(NULL != opt->key_seq) {
          tmap_error("options -1 and -2 cannot be used with -t", Exit, CommandLineArgument);
      }
      // OK
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
          switch(tmap_validate_flow_order(opt->flow_order)) {
            case 0:
              break;
            case -1:
              tmap_error("unrecognized DNA base (-x)", Exit, CommandLineArgument);
            case -2:
              tmap_error("all DNA bases must be present at least once (-x)", Exit, CommandLineArgument);
            default:
              tmap_error("unrecognized error (-x)", Exit, CommandLineArgument);
              break;
          }
      }
  }
  if(NULL != opt->key_seq) {
      if(0 == strcmp("sff", opt->key_seq) || 0 == strcmp("SFF", opt->key_seq)) {
          if(TMAP_READS_FORMAT_SFF != opt->reads_format) {
              tmap_error("an SFF was not specified (-r) but you want to use the sff flow order (-t)", Exit, CommandLineArgument);
          }
      }
      else {
          switch(tmap_validate_key_seq(opt->key_seq)) {
            case 0:
              break;
            case -1:
              tmap_error("unrecognized DNA base (-t)", Exit, CommandLineArgument); break;
              break;
            default:
              tmap_error("unrecognized error (-t)", Exit, CommandLineArgument);
              break;
          }
      }
  }
  tmap_error_cmd_check_int(opt->bw, 0, INT32_MAX, "-w");
  tmap_error_cmd_check_int(opt->softclip_type, 0, 3, "-g");
  if(1 == opt->softclip_key) {
      if(NULL == opt->key_seq && TMAP_READS_FORMAT_SFF != opt->reads_format) {
          tmap_error("an SFF (-r) or the key sequence (-t) must be specified to use -y", Exit, CommandLineArgument);
      }
  }
  else {
      tmap_error_cmd_check_int(opt->softclip_key, 0, 0, "-y");
  }

  tmap_error_cmd_check_int(opt->dup_window, -1, INT32_MAX, "-W");
  tmap_error_cmd_check_int(opt->max_seed_band, 1, INT32_MAX, "-B");
  tmap_error_cmd_check_int(opt->score_thr, INT32_MIN, INT32_MAX, "-T");
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
      tmap_error_cmd_check_int(opt->hit_frac, 0, 1, "-Y");
      break;
    case TMAP_MAP_ALGO_MAPALL:
      tmap_error_cmd_check_int(opt->aln_output_mode_ind, 0, 1, "-I");
      tmap_error_cmd_check_int(opt->mapall_score_thr, INT32_MIN, INT32_MAX, "-C");
      tmap_error_cmd_check_int(opt->mapall_mapq_thr, 0, 255, "-D");
      tmap_error_cmd_check_int(opt->mapall_keep_all, 0, 1, "-K");
      if(0 == opt->algos[0] || 0 == opt->num_stages) {
          tmap_error("no algorithms given for stage 1", Exit, CommandLineArgument);
      }
      for(i=0;i<2;i++) {
          // check mapping algorithm specific options
          tmap_map_opt_check(opt->opt_map1[i]);
          tmap_map_opt_check(opt->opt_map2[i]);
          tmap_map_opt_check(opt->opt_map3[i]);
          tmap_map_opt_check(opt->opt_map_vsw[i]);

          // check that common values match other opt values
          __tmap_map_opt_check_common1(opt, opt->opt_map1[i]);
          __tmap_map_opt_check_common1(opt, opt->opt_map2[i]);
          __tmap_map_opt_check_common1(opt, opt->opt_map3[i]);
          __tmap_map_opt_check_common1(opt, opt->opt_map_vsw[i]);
      }
      break;
    default:
      break;
  }
}

void
tmap_map_opt_print(tmap_map_opt_t *opt)
{
  int32_t i;
  fprintf(stderr, "algo_id=%d\n", opt->algo_id);
  fprintf(stderr, "fn_fasta=%s\n", opt->fn_fasta);
  for(i=0;i<opt->fn_reads_num;i++) {
      if(0 < i) fprintf(stderr, ",");
      fprintf(stderr, "fn_reads=%s\n", opt->fn_reads[i]);
  }
  fprintf(stderr, "fn_sam=%s\n", opt->fn_sam);
  fprintf(stderr, "reads_format=%d\n", opt->reads_format);
  fprintf(stderr, "reads_format=%d\n", opt->reads_format);
  fprintf(stderr, "score_match=%d\n", opt->score_match);
  fprintf(stderr, "pen_mm=%d\n", opt->pen_mm);
  fprintf(stderr, "pen_gapo=%d\n", opt->pen_gapo);
  fprintf(stderr, "pen_gape=%d\n", opt->pen_gape);
  fprintf(stderr, "fscore=%d\n", opt->fscore);
  fprintf(stderr, "flow_order=%s\n", opt->flow_order);
  fprintf(stderr, "flow_order_use_sff=%d\n", opt->flow_order_use_sff);
  fprintf(stderr, "key_seq=%s\n", opt->key_seq);
  fprintf(stderr, "key_seq_use_sff=%d\n", opt->key_seq_use_sff);
  fprintf(stderr, "bw=%d\n", opt->bw);
  fprintf(stderr, "softclip_type=%d\n", opt->softclip_type);
  fprintf(stderr, "softclip_key=%d\n", opt->softclip_key);
  fprintf(stderr, "dup_window=%d\n", opt->dup_window);
  fprintf(stderr, "max_seed_band=%d\n", opt->max_seed_band);
  fprintf(stderr, "score_thr=%d\n", opt->score_thr);
  fprintf(stderr, "reads_queue_size=%d\n", opt->reads_queue_size);
  fprintf(stderr, "num_threads=%d\n", opt->num_threads);
  fprintf(stderr, "aln_output_mode=%d\n", opt->aln_output_mode);
  fprintf(stderr, "sam_rg=%s\n", opt->sam_rg);
  fprintf(stderr, "sam_sff_tags=%d\n", opt->sam_sff_tags);
  fprintf(stderr, "remove_sff_clipping=%d\n", opt->remove_sff_clipping);
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
  fprintf(stderr, "hit_frac=%lf\n", opt->hit_frac);
  fprintf(stderr, "aln_output_mode_ind=%d\n", opt->aln_output_mode_ind);
  fprintf(stderr, "mapall_score_thr=%d\n", opt->mapall_score_thr);
  fprintf(stderr, "mapall_mapq_thr=%d\n", opt->mapall_mapq_thr);
  fprintf(stderr, "mapall_keep_all=%d\n", opt->mapall_keep_all);
}
