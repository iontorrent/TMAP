/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <config.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#include <unistd.h>
#endif
#include <unistd.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sam_print.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt_gen.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_bwt_match.h"
#include "../index/tmap_sa.h"
#include "../io/tmap_seq_io.h"
#include "../server/tmap_shm.h"
#include "../sw/tmap_sw.h"
#include "tmap_map_util.h"
#include "tmap_map1.h"
#include "tmap_map1_aux.h"
#include "tmap_map2.h"
#include "tmap_map2_aux.h"
#include "tmap_map3.h"
#include "tmap_map3_aux.h"
#include "tmap_map_driver.h"
#include "tmap_map_all.h"

/* Notes:
   We could avoid some computation give we know which algorithms to 
   run. This includes stacks, memory pools, as well as not loading 
   in all reference data and re-formatting all the input sequences.
   */

static void
tmap_map_all_sams_merge_helper(tmap_map_sams_t *dest, tmap_map_sams_t *src, int32_t stage)
{
  int32_t i, n;

  if(NULL == src || 0 == src->n) return;

  // make room
  n = dest->n;
  tmap_map_sams_realloc(dest, dest->n + src->n);
  // copy over
  for(i=0;i<src->n;i++) {
      // nullify
      tmap_map_sam_copy_and_nullify(&dest->sams[i+n], &src->sams[i]);
      // copy over the stage #
      dest->sams[i+n].algo_stage = stage;
  }
}

static tmap_map_sams_t *
tmap_map_all_sams_merge(tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2],
                        tmap_map_sams_t **sams_in, int32_t *algo_ids, int32_t n_algos,
                       int32_t stage, tmap_map_opt_t *opt)
{
  int32_t i;
  tmap_map_sams_t *sams = NULL;

  // init
  sams = tmap_map_sams_init();

  // remove duplicates before merging
  if(1 == opt->aln_output_mode_ind) {
      // smith waterman
      for(i=0;i<n_algos;i++) {
          tmap_map_util_sw(refseq, sams_in[i], seq, opt);
      }

      // duplicate removal
      for(i=0;i<n_algos;i++) {
          tmap_map_util_remove_duplicates(sams_in[i], opt->dup_window);
      }
  }

  // merge
  for(i=0;i<n_algos;i++) {
      tmap_map_all_sams_merge_helper(sams, sams_in[i], stage);
  }

  // no alignments
  if(0 == sams->n) return sams;

  // remove duplicates after merging
  if(0 == opt->aln_output_mode_ind) {
      // smith waterman
      sams = tmap_map_util_sw(refseq, sams, seq, opt);

      // duplicate removal
      tmap_map_util_remove_duplicates(sams, opt->dup_window);
  }

  // mapping quality
  tmap_map_util_mapq(sams, tmap_seq_get_bases(seq)->l, opt);

  // choose alignment(s)
  if(0 == opt->aln_output_mode_ind) {
      // consider all algos together
      tmap_map_sams_filter1(sams, opt->aln_output_mode, TMAP_MAP_ALGO_NONE);
  }
  else {
      // consider all algos independently
      for(i=0;i<n_algos;i++) {
          tmap_map_sams_filter1(sams, opt->aln_output_mode, algo_ids[i]);
      }
  }

  return sams;
}

// thread data
typedef struct {
    void *data_map1[2];
    void *data_map2[2];
    void *data_map3[2];
} tmap_map_all_thread_data_t;

static int32_t
tmap_map_all_init(tmap_refseq_t *refseq, tmap_map_opt_t *opt)
{
  int32_t i;

  for(i=0;i<opt->num_stages;i++) {
      // map1
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP1) {
          opt->opt_map1[i]->algo_stage = i+1;
          tmap_map1_init(refseq, opt->opt_map1[i]);
      }
      // map2
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP2) {
          opt->opt_map2[i]->algo_stage = i+1;
          tmap_map2_init(refseq, opt->opt_map2[i]);
      }
      // map3
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP3) {
          opt->opt_map3[i]->algo_stage = i+1;
          tmap_map3_init(refseq, opt->opt_map3[i]);
      }
  }
  return 0;
}

static int32_t
tmap_map_all_thread_init(void **data, tmap_map_opt_t *opt)
{
  int32_t i;
  tmap_map_all_thread_data_t *d = NULL;

  d = tmap_calloc(1, sizeof(tmap_map_all_thread_data_t), "d");
  
  for(i=0;i<opt->num_stages;i++) {
      // map1
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP1) {
          tmap_map1_thread_init(&d->data_map1[i], opt->opt_map1[i]);
      }
      // map2
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP2) {
          tmap_map2_thread_init(&d->data_map2[i], opt->opt_map2[i]);
      }
      // map3
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP3) {
          tmap_map3_thread_init(&d->data_map3[i], opt->opt_map3[i]);
      }
  }

  (*data) = (void*)d;

  return 0;
}

// Note: we assume that filtering and duplicate removal will not be done but the
// map* algorithms when their program options "algo_id != -1"
static tmap_map_sams_t*
tmap_map_all_thread_map(void **data, tmap_seq_t *seq, tmap_refseq_t *refseq, tmap_bwt_t *bwt[2], tmap_sa_t *sa[2], tmap_map_opt_t *opt)
{
  int32_t i, n_algos=3;
  tmap_map_sams_t *sams = NULL;
  tmap_map_sams_t *sams_in[3] = {NULL,NULL,NULL};
  int32_t algo_ids[3] = {TMAP_MAP_ALGO_MAP1,TMAP_MAP_ALGO_MAP2,TMAP_MAP_ALGO_MAP3};
  tmap_map_all_thread_data_t *d = (tmap_map_all_thread_data_t*)(*data);
  
  // for mapall optimization
  int32_t seq_len;
  tmap_seq_t *seqs[4]={NULL,NULL,NULL,NULL};
  tmap_seq_t *seqs_tmp[2]={NULL,NULL};
  tmap_string_t *bases[4]={NULL,NULL,NULL,NULL};
  tmap_string_t *bases_tmp[2]={NULL,NULL};

  seq_len = tmap_seq_get_bases(seq)->l; 

  for(i=0;i<4;i++) {
      // TODO: only if necessary
      seqs[i]= tmap_seq_clone(seq); // clone the sequence 
      switch(i) {
        case 0: // forward
          break;
        case 1: // reverse compliment
          tmap_seq_reverse_compliment(seqs[i]); break;
        case 2: // reverse
          tmap_seq_reverse(seqs[i]); break;
        case 3: // compliment
          tmap_seq_compliment(seqs[i]); break;
      }
      tmap_seq_to_int(seqs[i]); // convert to integers
      bases[i] = tmap_seq_get_bases(seqs[i]); // get bases
  }

  for(i=0;i<opt->num_stages;i++) {
      // map1
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP1) {
          seqs_tmp[0] = seqs[2]; seqs_tmp[1] = seqs[1];
          bases_tmp[0] = bases[2]; bases_tmp[1] = bases[1];
          sams_in[0] = tmap_map1_thread_map_core(&d->data_map1[i], seqs_tmp, bases_tmp, seq_len, refseq, bwt, sa, opt->opt_map1[i]);
      }
      else {
          sams_in[0] = tmap_map_sams_init();
      }
      // map2
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP2) {
          sams_in[1] = tmap_map2_thread_map_core(&d->data_map2[i], seqs, seq_len, refseq, bwt, sa, opt->opt_map2[i]);
      }
      else {
          sams_in[1] = tmap_map_sams_init();
      }
      // map3
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP3) {
          seqs_tmp[0] = seqs[0];
          seqs_tmp[1] = seqs[1];
          sams_in[2] = tmap_map3_thread_map_core(&d->data_map3[i], seqs, seq_len, refseq, bwt, sa, opt->opt_map3[i]);
      }
      else {
          sams_in[2] = tmap_map_sams_init();
      }

      // merge all the mappings
      sams = tmap_map_all_sams_merge(seq, refseq, bwt, sa,
                                     sams_in, algo_ids, n_algos, 
                                     i+1, opt);
      for(i=0;i<n_algos;i++) {
          tmap_map_sams_destroy(sams_in[i]);
          sams_in[i] = NULL;
      }

      // did we find any mappings?
      if(0 < sams->n) {
          break;
      }
      else if(0 == i && 2 == opt->num_stages) {
          // destroy for the next round
          tmap_map_sams_destroy(sams);
          sams = NULL;
      }
  }
  
  // destroy
  for(i=0;i<4;i++) {
      tmap_seq_destroy(seqs[i]);
  }

  return sams;
}

static int32_t
tmap_map_all_thread_cleanup(void **data, tmap_map_opt_t *opt)
{
  int32_t i;
  tmap_map_all_thread_data_t *d = (tmap_map_all_thread_data_t*)(*data);

  for(i=0;i<opt->num_stages;i++) {
      // map1
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP1) {
          tmap_map1_thread_cleanup(&d->data_map1[i], opt->opt_map1[i]);
      }
      // map2
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP2) {
          tmap_map2_thread_cleanup(&d->data_map2[i], opt->opt_map2[i]);
      }
      // map3
      if(opt->algos[i] & TMAP_MAP_ALGO_MAP3) {
          tmap_map3_thread_cleanup(&d->data_map3[i], opt->opt_map3[i]);
      }
  }
  free(d);
  (*data) = NULL;
  return 0;
}

static void
tmap_map_all_core(tmap_map_opt_t *opt)
{
  // run the driver
  tmap_map_driver_core(tmap_map_all_init,
                       tmap_map_all_thread_init,
                       tmap_map_all_thread_map,
                       NULL,
                       tmap_map_all_thread_cleanup,
                       opt);
}

// for map1/map2/map3
#define __tmap_map_all_opts_copy1(opt_map_all, opt_map_other) do { \
    (opt_map_other)->fn_fasta = tmap_strdup((opt_map_all)->fn_fasta); \
    (opt_map_other)->fn_reads = tmap_strdup((opt_map_all)->fn_reads); \
    (opt_map_other)->reads_format = (opt_map_all)->reads_format; \
    (opt_map_other)->score_match = (opt_map_all)->score_match; \
    (opt_map_other)->pen_mm = (opt_map_all)->pen_mm; \
    (opt_map_other)->pen_gapo = (opt_map_all)->pen_gapo; \
    (opt_map_other)->pen_gape = (opt_map_all)->pen_gape; \
    (opt_map_other)->fscore = (opt_map_all)->fscore; \
    (opt_map_other)->flow_order = tmap_strdup((opt_map_all)->flow_order); \
    (opt_map_other)->flow_order_use_sff = (opt_map_all)->flow_order_use_sff; \
    (opt_map_other)->bw = (opt_map_all)->bw; \
    (opt_map_other)->softclip_type = (opt_map_all)->softclip_type; \
    (opt_map_other)->dup_window = (opt_map_all)->dup_window; \
    (opt_map_other)->max_seed_band = (opt_map_all)->max_seed_band; \
    (opt_map_other)->score_thr = (opt_map_all)->score_thr; \
    (opt_map_other)->reads_queue_size = (opt_map_all)->reads_queue_size; \
    (opt_map_other)->num_threads = (opt_map_all)->num_threads; \
    (opt_map_other)->aln_output_mode = TMAP_MAP_UTIL_ALN_MODE_ALL; \
    (opt_map_other)->sam_rg = tmap_strdup((opt_map_all)->sam_rg); \
    (opt_map_other)->sam_sff_tags = (opt_map_all)->sam_sff_tags; \
    (opt_map_other)->input_compr = (opt_map_all)->input_compr; \
    (opt_map_other)->output_compr = (opt_map_all)->output_compr; \
    (opt_map_other)->shm_key = (opt_map_all)->shm_key; \
} while(0)

int32_t
tmap_map_all_opt_parse(int argc, char *argv[], tmap_map_opt_t *opt)
{
  int32_t i, j, start, opt_type, opt_type_next, opt_stage, opt_stage_next;

  // parse common options as well as map1/map2/map3 commands
  start = 0;
  i = 1;
  opt_type = opt_type_next = TMAP_MAP_ALGO_NONE;
  opt_stage = opt_stage_next = 0;
  opt->num_stages = 0;
  while(i<argc) {
      if(0 == strcmp("map1", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP1;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
          if(0 == opt->num_stages) opt->num_stages++;
      }
      else if(0 == strcmp("map2", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP2;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
          if(0 == opt->num_stages) opt->num_stages++;
      }
      else if(0 == strcmp("map3", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP3;
          opt->algos[0] |= opt_type_next;
          opt_stage_next = 1;
          if(0 == opt->num_stages) opt->num_stages++;
      }
      else if(0 == strcmp("MAP1", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP1;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
          if(1 == opt->num_stages) opt->num_stages++;
      }
      else if(0 == strcmp("MAP2", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP2;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
          if(1 == opt->num_stages) opt->num_stages++;
      }
      else if(0 == strcmp("MAP3", argv[i])) {
          opt_type_next = TMAP_MAP_ALGO_MAP3;
          opt->algos[1] |= opt_type_next;
          opt_stage_next = 2;
          if(1 == opt->num_stages) opt->num_stages++;
      }

      /*
         fprintf(stderr, "i=%d start=%d argc=%d opt_type=%d opt_type_next=%d argv[%d]=%s\n",
         i, start, argc, opt_type, opt_type_next, i, argv[i]);
         */

      if(opt_type != opt_type_next
         || i == argc-1) {
          if(i == argc-1) {
              i++;
          }
          optind=1; // needed for getopt
          switch(opt_type) {
            case TMAP_MAP_ALGO_NONE:
              // parse common options
              if(0 == tmap_map_opt_parse(i-start, argv+start, opt)) {
                  return 0;
              }
              // copy over common values into the other opts
              for(j=0;j<2;j++) {
                  __tmap_map_all_opts_copy1(opt, opt->opt_map1[j]);
                  __tmap_map_all_opts_copy1(opt, opt->opt_map2[j]);
                  __tmap_map_all_opts_copy1(opt, opt->opt_map3[j]);
              }
              break;
            case TMAP_MAP_ALGO_MAP1:
              if(0 < i - start) {
                  if(0 == tmap_map_opt_parse(i-start, argv+start, opt->opt_map1[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            case TMAP_MAP_ALGO_MAP2:
              // parse map2 options
              if(0 < i - start) {
                  if(0 == tmap_map_opt_parse(i-start, argv+start, opt->opt_map2[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            case TMAP_MAP_ALGO_MAP3:
              // parse map3 options
              if(0 < i - start) {
                  if(0 == tmap_map_opt_parse(i-start, argv+start, opt->opt_map3[opt_stage-1])) {
                      return 0;
                  }
              }
              break;
            default:
              tmap_error("bug encountered", Exit, OutOfRange);
          }
          opt_type = opt_type_next;
          opt_stage = opt_stage_next;
          start = i;
      }
      i++;
      if(argc < i) {
          i = argc;
      }
  }
  optind = i;
  
  // do this after parsing
  opt->argc = argc; opt->argv = argv;

  return 1;
}

int 
tmap_map_all_main(int argc, char *argv[])
{
  tmap_map_opt_t *opt = NULL;

  // init opt
  opt = tmap_map_opt_init(TMAP_MAP_ALGO_MAPALL);
      
  // get options
  if(1 != tmap_map_all_opt_parse(argc, argv, opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return tmap_map_opt_usage(opt);
  }
  else { 
      // check command line arguments
      tmap_map_opt_check(opt);
  }

  // run map_all
  tmap_map_all_core(opt);

  // destroy opt
  tmap_map_opt_destroy(opt);

  tmap_progress_print2("terminating successfully");

  return 0;
}
