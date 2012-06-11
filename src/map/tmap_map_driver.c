/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <config.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif
#include <unistd.h>
#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_definitions.h"
#include "../util/tmap_progress.h"
#include "../util/tmap_sam_convert.h"
#include "../util/tmap_sort.h"
#include "../util/tmap_rand.h"
#include "../util/tmap_hash.h"
#include "../seq/tmap_seq.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt_gen.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_bwt_match.h"
#include "../index/tmap_bwt_match_hash.h"
#include "../index/tmap_sa.h"
#include "../index/tmap_index.h"
#include "../io/tmap_seqs_io.h"
#include "../server/tmap_shm.h"
#include "../sw/tmap_fsw.h"
#include "../sw/tmap_sw.h"
#include "util/tmap_map_stats.h"
#include "util/tmap_map_util.h"
#include "pairing/tmap_map_pairing.h"
#include "tmap_map_driver.h"

// NB: do not turn these on, as they do not currently improve run time. They
// could be useful if many duplicate lookups are performed and the hash
// retrieval was fast...
//#define TMAP_DRIVER_USE_HASH 1
//#define TMAP_DRIVER_CLEAR_HASH_PER_READ 1

#define __tmap_map_sam_sort_score_lt(a, b) ((a).score > (b).score)
TMAP_SORT_INIT(tmap_map_sam_sort_score, tmap_map_sam_t, __tmap_map_sam_sort_score_lt)

static void
tmap_map_driver_do_init(tmap_map_driver_t *driver, tmap_refseq_t *refseq)
{
  int32_t i, j;
  for(i=0;i<driver->num_stages;i++) {
      tmap_map_driver_stage_t *stage = driver->stages[i];
      for(j=0;j<stage->num_algorithms;j++) {
          tmap_map_driver_algorithm_t *algorithm = stage->algorithms[j];
          if(NULL != algorithm->func_init 
             && 0 != algorithm->func_init(&algorithm->data, refseq, algorithm->opt)) {
              tmap_error("the thread function could not be initialized", Exit, OutOfRange);
          }
      }
  }
}

static void
tmap_map_driver_do_cleanup(tmap_map_driver_t *driver)
{
  int32_t i, j;
  for(i=0;i<driver->num_stages;i++) {
      tmap_map_driver_stage_t *stage = driver->stages[i];
      for(j=0;j<stage->num_algorithms;j++) {
          tmap_map_driver_algorithm_t *algorithm = stage->algorithms[j];
          if(NULL != algorithm->func_cleanup 
             && 0 != algorithm->func_cleanup(&algorithm->data)) {
              tmap_error("the thread function could not be cleaned up", Exit, OutOfRange);
          }
      }
  }
}

static void
tmap_map_driver_do_threads_init(tmap_map_driver_t *driver, 
                                int32_t tid)
{
  int32_t i, j;
  for(i=0;i<driver->num_stages;i++) {
      tmap_map_driver_stage_t *stage = driver->stages[i];
      for(j=0;j<stage->num_algorithms;j++) {
          tmap_map_driver_algorithm_t *algorithm = stage->algorithms[j];
          if(NULL != algorithm->func_thread_init 
             && 0 != algorithm->func_thread_init(&algorithm->thread_data[tid], 
                                                 algorithm->opt)) {
              tmap_error("the thread function could not be initialized", Exit, OutOfRange);
          }
      }
  }
}

static void
tmap_map_driver_do_threads_cleanup(tmap_map_driver_t *driver, int32_t tid)
{
  int32_t i, j;
  for(i=0;i<driver->num_stages;i++) {
      tmap_map_driver_stage_t *stage = driver->stages[i];
      for(j=0;j<stage->num_algorithms;j++) {
          tmap_map_driver_algorithm_t *algorithm = stage->algorithms[j];
          if(NULL != algorithm->func_thread_cleanup 
             && 0 != algorithm->func_thread_cleanup(&algorithm->thread_data[tid], algorithm->opt)) {
              tmap_error("the thread function could not be cleaned up", Exit, OutOfRange);
          }
      }
  }
}

static void
tmap_map_driver_init_seqs(tmap_seq_t **seqs, tmap_seq_t *seq, int32_t max_length)
{
  int32_t i;
  // init seqs
  for(i=0;i<4;i++) {
      // TODO: only if necessary
      seqs[i] = tmap_seq_clone(seq); // clone the sequence 
      // modify the length before reversing or reverse complimenting 
      if(0 < max_length && max_length < tmap_seq_get_bases_length(seq)) { // NB: does not modify quality string or other meta data
          tmap_seq_get_bases(seqs[i])->l = max_length;
          tmap_seq_get_bases(seqs[i])->s[max_length] = '\0';
      }
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
  }
}

void
tmap_map_driver_core_worker(sam_header_t *sam_header,
                            tmap_seqs_t **seqs_buffer, 
                            tmap_map_record_t **records, 
                            tmap_map_bams_t **bams,
                            int32_t seqs_buffer_length,
                            int32_t *buffer_idx,
                            tmap_index_t *index,
                            tmap_map_driver_t *driver,
                            tmap_map_stats_t *stat,
                            tmap_rand_t *rand,
                            int32_t tid)
{
  int32_t i, j, k, low = 0;
  int32_t found;
  tmap_seq_t ***seqs = NULL;
  tmap_bwt_match_hash_t *hash=NULL;
  int32_t max_num_ends = 0;

#ifdef TMAP_DRIVER_USE_HASH
  // init the occurence hash
  hash = tmap_bwt_match_hash_init(); 
#endif

  // init memory
  max_num_ends = 2;
  seqs = tmap_malloc(sizeof(tmap_seq_t**)*max_num_ends, "seqs");
  for(i=0;i<max_num_ends;i++) {
      seqs[i] = tmap_calloc(4, sizeof(tmap_seq_t*), "seqs[i]");
  }
  
  // initialize thread data
  tmap_map_driver_do_threads_init(driver, tid);

  // Go through the buffer
  while(low < seqs_buffer_length) {
      if(tid == (low % driver->opt->num_threads)) {
          tmap_map_stats_t *stage_stat = NULL;
          tmap_map_record_t *record_prev = NULL;
          int32_t num_ends;
          
#ifdef TMAP_DRIVER_USE_HASH
#ifdef TMAP_DRIVER_CLEAR_HASH_PER_READ
          // TODO: should we hash each read, or across the thread?
          tmap_bwt_match_hash_clear(hash);
#endif
#endif

          num_ends = seqs_buffer[low]->n;
          if(max_num_ends < num_ends) {
              seqs = tmap_realloc(seqs, sizeof(tmap_seq_t**)*num_ends, "seqs");
              while(max_num_ends < num_ends) {
                  seqs[max_num_ends] = tmap_calloc(4, sizeof(tmap_seq_t*), "seqs[max_num_ends]");
                  max_num_ends++;
              }
              max_num_ends = num_ends;
          }
          
          // re-initialize the random seed
          if(driver->opt->rand_read_name) {
              tmap_rand_reinit(rand, tmap_hash_str_hash_func(tmap_seq_get_name(seqs_buffer[low]->seqs[0])->s));
          }

          // init
          for(i=0;i<num_ends;i++) {
              tmap_map_driver_init_seqs(seqs[i], seqs_buffer[low]->seqs[i], -1);
              stat->num_reads++;
          }

          // init records
          records[low] = tmap_map_record_init(num_ends);

          // go through each stage
          for(i=0;i<driver->num_stages;i++) { // for each stage
              tmap_map_driver_stage_t *stage = driver->stages[i];

              // stage stats
              stage_stat = tmap_map_stats_init();

              // seed
              for(j=0;j<num_ends;j++) { // for each end
                  tmap_seq_t **stage_seqs = NULL;
                  // should we seed using the whole read?
                  if(0 < stage->opt->stage_seed_max_length && stage->opt->stage_seed_max_length < tmap_seq_get_bases_length(seqs[j][0])) {
                      stage_seqs = tmap_calloc(4, sizeof(tmap_seq_t*), "seqs[i]");
                      tmap_map_driver_init_seqs(stage_seqs, seqs_buffer[low]->seqs[i], stage->opt->stage_seed_max_length);
                  }
                  else {
                      stage_seqs = seqs[j];
                  }
                  for(k=0;k<stage->num_algorithms;k++) { // for each algorithm
                      tmap_map_driver_algorithm_t *algorithm = stage->algorithms[k];
                      tmap_map_sams_t *sams = NULL;
                      if(i+1 != algorithm->opt->algo_stage) {
                          tmap_bug();
                      }
                      // map
                      sams = algorithm->func_thread_map(&algorithm->thread_data[tid], stage_seqs, index, hash, rand, algorithm->opt);
                      if(NULL == sams) {
                          tmap_error("the thread function did not return a mapping", Exit, OutOfRange);
                      }
                      // append
                      tmap_map_sams_merge(records[low]->sams[j], sams);
                      // destroy
                      tmap_map_sams_destroy(sams);
                  }
                  stage_stat->num_after_seeding += records[low]->sams[j]->n;
                  if(0 < stage->opt->stage_seed_max_length && stage->opt->stage_seed_max_length < tmap_seq_get_bases_length(seqs[j][0])) {
                      // free
                      for(j=0;j<4;j++) {
                          tmap_seq_destroy(stage_seqs[j]);
                          stage_seqs[j] = NULL;
                      }
                  }
                  stage_seqs = NULL; // do not use
              }

              // keep mappings for subsequent stages or restore mappings from
              // previous stages
              if(1 == stage->opt->stage_keep_all) {
                  // merge from the previous stage
                  if(0 < i) {
                      tmap_map_record_merge(records[low], record_prev);
                      // destroy the record
                      tmap_map_record_destroy(record_prev);
                      record_prev = NULL;
                  }

                  // keep for the next stage
                  if(i < driver->num_stages-1) { // more stages left
                      record_prev = tmap_map_record_clone(records[low]);
                  }
              }

              // generate scores with smith waterman
              for(j=0;j<num_ends;j++) { // for each end
                  records[low]->sams[j] = tmap_map_util_sw_gen_score(index->refseq, records[low]->sams[j], seqs[j], rand, stage->opt);
                  stage_stat->num_after_scoring += records[low]->sams[j]->n;
              }

              // remove duplicates
              for(j=0;j<num_ends;j++) { // for each end
                  tmap_map_util_remove_duplicates(records[low]->sams[j], stage->opt->dup_window, rand);
                  stage_stat->num_after_rmdup += records[low]->sams[j]->n;
              }
              
              // (single-end) mapping quality
              for(j=0;j<num_ends;j++) { // for each end
                  driver->func_mapq(records[low]->sams[j], tmap_seq_get_bases_length(seqs[j][0]), stage->opt);
              }

              // filter if we have more stages
              if(i < driver->num_stages-1) {
                  for(j=0;j<num_ends;j++) { // for each end
                      tmap_map_sams_filter2(records[low]->sams[j], stage->opt->stage_score_thr, stage->opt->stage_mapq_thr);
                  }
              }

              if(0 <= driver->opt->strandedness && 0 <= driver->opt->positioning
                 && 2 == num_ends && 0 < records[low]->sams[0]->n && 0 < records[low]->sams[1]->n) { // pairs of reads!

                  // read rescue
                  if(1 == stage->opt->read_rescue) {
                      int32_t flag = tmap_map_pairing_read_rescue(index->refseq, 
                                                                  records[low]->sams[0], records[low]->sams[1],
                                                                  seqs[0], seqs[1],
                                                                  rand, stage->opt);
                      // recalculate mapping qualities if necessary
                      if(0 < (flag & 0x1)) { // first end was rescued
                          //fprintf(stderr, "re-doing mapq for end #1\n");
                          driver->func_mapq(records[low]->sams[0], tmap_seq_get_bases_length(seqs[0][0]), stage->opt);
                      }
                      if(0 < (flag & 0x2)) { // second end was rescued
                          //fprintf(stderr, "re-doing mapq for end #2\n");
                          driver->func_mapq(records[low]->sams[1], tmap_seq_get_bases_length(seqs[1][0]), stage->opt);
                      }
                  }
                  // pick pairs
                  tmap_map_pairing_pick_pairs(records[low]->sams[0], records[low]->sams[1],
                                              seqs[0][0], seqs[1][0],
                                              rand, stage->opt);
                  // TODO: if we have one end for a pair, do we go onto the second
                  // stage?
              }
              else {

                  // choose alignments
                  for(j=0;j<num_ends;j++) { // for each end
                      tmap_map_sams_filter1(records[low]->sams[j], stage->opt->aln_output_mode, TMAP_MAP_ALGO_NONE, rand);
                      stage_stat->num_after_filter += records[low]->sams[j]->n;
                  }
              }

              // generate the cigars
              found = 0;
              for(j=0;j<num_ends;j++) { // for each end
                  records[low]->sams[j] = tmap_map_util_sw_gen_cigar(index->refseq, records[low]->sams[j], seqs[j], stage->opt);
                  if(0 < records[low]->sams[j]->n) {
                      stage_stat->num_with_mapping++;
                      found = 1;
                  }
              }

              // TODO
              // if paired, update pairing score based on target start?

              // did we find any mappings?
              if(1 == found) { // yes
                  tmap_map_stats_add(stat, stage_stat);
                  tmap_map_stats_destroy(stage_stat);
                  break;
              }
              else { // no
                  tmap_map_record_destroy(records[low]);
                  // re-init
                  records[low] = tmap_map_record_init(num_ends);
              }
              tmap_map_stats_destroy(stage_stat);
          }

          // flowspace re-align and sorting
          if(1 == driver->opt->aln_flowspace) {
              for(i=0;i<num_ends;i++) {
                  if(0 < records[low]->sams[i]->n) {
                      //stat->num_with_mapping++;
                      // re-align the alignments in flow-space
                      tmap_seq_t *seq = seqs_buffer[low]->seqs[i];
                      // TODO: if this is run, we do not need to run tmap_sw_global_banded_core...
                      // NB: seqs_buffer should have its key sequence if 0 < key_seq_len
                      if(1 == tmap_map_util_fsw(seq, records[low]->sams[i], index->refseq, 
                                        driver->opt->bw, driver->opt->softclip_type, driver->opt->score_thr,
                                        driver->opt->score_match, driver->opt->pen_mm, driver->opt->pen_gapo,
                                        driver->opt->pen_gape, driver->opt->fscore, 1-driver->opt->ignore_flowgram)) {
                          // sort by alignment score
                          if(1 < records[low]->sams[i]->n) {
                              tmap_sort_introsort(tmap_map_sam_sort_score,
                                                  records[low]->sams[i]->n, 
                                                  records[low]->sams[i]->sams);
                          }
                      }
                  }
                  if(NULL == seqs_buffer[low]->seqs[i]) {
                      tmap_error("bug encoutereed", Exit, OutOfRange);
                  }
              }
          }

          // convert the record to bam
          if(1 == seqs_buffer[low]->n) {
              bams[low] = tmap_map_bams_init(1); 
              bams[low]->bams[0] = tmap_map_sams_print(seqs_buffer[low]->seqs[0], index->refseq, records[low]->sams[0], 
                                                          0, NULL, driver->opt->sam_flowspace_tags, driver->opt->bidirectional, driver->opt->seq_eq);
          }
          else {
              bams[low] = tmap_map_bams_init(seqs_buffer[low]->n);
              for(j=0;j<seqs_buffer[low]->n;j++) {
                  bams[low]->bams[j] = tmap_map_sams_print(seqs_buffer[low]->seqs[j], index->refseq, records[low]->sams[j],
                                                           (0 == j) ? 1 : ((seqs_buffer[low]->n-1 == j) ? 2 : 0),
                                                           records[low]->sams[(j+1) % seqs_buffer[low]->n], 
                                                           driver->opt->sam_flowspace_tags, driver->opt->bidirectional, driver->opt->seq_eq);
              }
          }

          // free alignments, for space
          tmap_map_record_destroy(records[low]); 
          records[low] = NULL;

          // free seqs
          for(i=0;i<num_ends;i++) {
              for(j=0;j<4;j++) {
                  tmap_seq_destroy(seqs[i][j]);
                  seqs[i][j] = NULL;
              }
          }
          tmap_map_record_destroy(record_prev);
      }
      // next
      (*buffer_idx) = low;
      low++;
  }
  (*buffer_idx) = seqs_buffer_length;
                  
  // free thread variables
  for(i=0;i<max_num_ends;i++) {
      free(seqs[i]);
  }
  free(seqs);

  // cleanup
  tmap_map_driver_do_threads_cleanup(driver, tid);
#ifdef TMAP_DRIVER_USE_HASH
  // free hash
  tmap_bwt_match_hash_destroy(hash);
#endif
}

void *
tmap_map_driver_core_thread_worker(void *arg)
{
  tmap_map_driver_thread_data_t *thread_data = (tmap_map_driver_thread_data_t*)arg;

  tmap_map_driver_core_worker(thread_data->sam_header, thread_data->seqs_buffer, thread_data->records, thread_data->bams, 
                              thread_data->seqs_buffer_length, thread_data->buffer_idx, thread_data->index, thread_data->driver, 
                              thread_data->stat, thread_data->rand, thread_data->tid);

  return arg;
}

void 
tmap_map_driver_core(tmap_map_driver_t *driver)
{
  uint32_t i, j, k, n_reads_processed=0; // # of reads processed
  int32_t seqs_buffer_length=0; // # of reads read in
  tmap_seqs_io_t *io_in = NULL; // input file(s)
  tmap_sam_io_t *io_out = NULL; // output file
  tmap_seqs_t **seqs_buffer = NULL; // buffer for the reads
  tmap_map_record_t **records=NULL; // buffer for the mapped data
  tmap_map_bams_t **bams=NULL;// buffer for the mapped BAM data
  tmap_index_t *index = NULL; // reference indes
  tmap_map_stats_t *stat = NULL; // alignment statistics
#ifdef HAVE_LIBPTHREAD
  pthread_attr_t *attr = NULL;
  pthread_t *threads = NULL;
  tmap_map_driver_thread_data_t *thread_data=NULL;
  tmap_rand_t **rand = NULL; // random # generator for each thread
  tmap_map_stats_t **stats = NULL; // alignment statistics for each thread
#else
  tmap_rand_t *rand = NULL; // random # generator
#endif
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
  tmap_rand_t *rand_core = tmap_rand_init(13); // random # generator for sampling
  uint32_t k;
#endif
  int32_t seq_type, reads_queue_size; // read type, read queue size
  bam_header_t *header = NULL; // BAM Header
  int32_t buffer_idx; // buffer index for processing data with a single thread


  /*
  if(NULL == driver->opt->fn_reads) {
      tmap_progress_set_verbosity(0); 
  }
  */
          
  tmap_progress_print("running with %d threads (%s)",
                       driver->opt->num_threads,
                       (0 == driver->opt->num_threads_autodetected) ? "user set" : "autodetected");
  
  // print out the algorithms and stages
  for(i=0;i<driver->num_stages;i++) {
      for(j=0;j<driver->stages[i]->num_algorithms;j++) {
          tmap_progress_print("%s will be run in stage %d", 
                               tmap_algo_id_to_name(driver->stages[i]->algorithms[j]->opt->algo_id),
                               driver->stages[i]->algorithms[j]->opt->algo_stage);
      }
  }

  // open the reads file for reading
  // NB: may have no fns (streaming in)
  seq_type = tmap_reads_format_to_seq_type(driver->opt->reads_format); 
  io_in = tmap_seqs_io_init(driver->opt->fn_reads, driver->opt->fn_reads_num, seq_type, driver->opt->input_compr);

  // get the index
  index = tmap_index_init(driver->opt->fn_fasta, driver->opt->shm_key);

  // initialize the driver->options and print any relevant information
  tmap_map_driver_do_init(driver, index->refseq);

  // allocate the buffer
  if(-1 == driver->opt->reads_queue_size) {
      reads_queue_size = 1;
  }
  else {
      reads_queue_size = driver->opt->reads_queue_size;
  }
  seqs_buffer = tmap_malloc(sizeof(tmap_seqs_t*)*reads_queue_size, "seqs_buffer");
  for(i=0;i<reads_queue_size;i++) { // initialize the buffer
      seqs_buffer[i] = tmap_seqs_init(seq_type);
  }
  records = tmap_malloc(sizeof(tmap_map_record_t*)*reads_queue_size, "records");
  bams = tmap_malloc(sizeof(tmap_map_bams_t*)*reads_queue_size, "bams");

  stat = tmap_map_stats_init();
#ifdef HAVE_LIBPTHREAD
  stats = tmap_malloc(driver->opt->num_threads * sizeof(tmap_map_stats_t*), "stats");
  rand = tmap_malloc(driver->opt->num_threads * sizeof(tmap_rand_t*), "rand");
  for(i=0;i<driver->opt->num_threads;i++) {
      stats[i] = tmap_map_stats_init();
      rand[i] = tmap_rand_init(i);
  }
#else
  rand = tmap_rand_init(13);
#endif
  
  // BAM Header
  header = tmap_seqs_io_to_bam_header(index->refseq, io_in, 
                                      driver->opt->sam_rg, driver->opt->sam_rg_num,
                                      driver->opt->argc, driver->opt->argv);

  // open the output file
  switch(driver->opt->output_type) {
    case 0: // SAM
      io_out = tmap_sam_io_init2((NULL == driver->opt->fn_sam) ? "-" : driver->opt->fn_sam, "wh", header); 
      break;
    case 1:
      io_out = tmap_sam_io_init2((NULL == driver->opt->fn_sam) ? "-" : driver->opt->fn_sam, "wb", header); 
      break;
    case 2:
      io_out = tmap_sam_io_init2((NULL == driver->opt->fn_sam) ? "-" : driver->opt->fn_sam, "wbu", header); 
      break;
    default:
      tmap_bug();
  }

  // destroy the BAM Header
  bam_header_destroy(header);
  header = NULL;

  // main processing loop
  tmap_progress_print("processing reads");
  while(1) {
      tmap_progress_print("loading reads");
      // get the reads
      seqs_buffer_length = tmap_seqs_io_read_buffer(io_in, seqs_buffer, reads_queue_size, io_out->fp->header->header);
      tmap_progress_print2("loaded %d reads", seqs_buffer_length);
      if(0 == seqs_buffer_length) {
          break;
      }
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
      // sample reads
      if(driver->opt->sample_reads < 1) {
          for(i=j=0;i<seqs_buffer_length;i++) {
              if(driver->opt->sample_reads < tmap_rand_get(rand_core)) continue; // skip
              if(j < i) {
                  tmap_seqs_t *seqs;
                  seqs = seqs_buffer[j];
                  seqs_buffer[j] = seqs_buffer[i]; 
                  seqs_buffer[i] = seqs;
              }
              j++;
          }
          tmap_progress_print2("sampling %d out of %d [%.2lf%%]", j, seqs_buffer_length, 100.0*j/(double)seqs_buffer_length);
          seqs_buffer_length = j;
          if(0 == seqs_buffer_length) continue;
      }
#endif

      // do alignment
#ifdef HAVE_LIBPTHREAD
      if(1 == driver->opt->num_threads) {
          buffer_idx = 0;
          tmap_map_driver_core_worker(io_out->fp->header->header, seqs_buffer, records, bams, 
                                      seqs_buffer_length, &buffer_idx, index, driver, stat, rand[0], 0);
      }
      else {
          attr = tmap_calloc(1, sizeof(pthread_attr_t), "attr");
          pthread_attr_init(attr);
          pthread_attr_setdetachstate(attr, PTHREAD_CREATE_JOINABLE);

          threads = tmap_calloc(driver->opt->num_threads, sizeof(pthread_t), "threads");
          thread_data = tmap_calloc(driver->opt->num_threads, sizeof(tmap_map_driver_thread_data_t), "thread_data");

          // create threads
          for(i=0;i<driver->opt->num_threads;i++) {
              thread_data[i].sam_header = io_out->fp->header->header;
              thread_data[i].seqs_buffer = seqs_buffer;
              thread_data[i].seqs_buffer_length = seqs_buffer_length;
              thread_data[i].buffer_idx = tmap_calloc(1, sizeof(int32_t), "thread_data[i].buffer_id");
              thread_data[i].records = records;
              thread_data[i].bams = bams;
              thread_data[i].index = index;
              thread_data[i].driver = driver;
              thread_data[i].stat = stats[i];
              thread_data[i].rand = rand[i];
              thread_data[i].tid = i;
              if(0 != pthread_create(&threads[i], attr, tmap_map_driver_core_thread_worker, &thread_data[i])) {
                  tmap_error("error creating threads", Exit, ThreadError);
              }
          }
      }
#else 
      buffer_idx = 0;
      tmap_map_driver_core_worker(io_out->fp->header->header, seqs_buffer, records, bams, 
                                  seqs_buffer_length, &buffer_idx, index, driver, stat, rand, 0);
#endif

      /*
      if(-1 != driver->opt->reads_queue_size) {
          tmap_progress_print("writing alignments");
      }
      */

      // write data
      for(i=0;i<seqs_buffer_length;i++) {
#ifdef HAVE_LIBPTHREAD
          if(1 < driver->opt->num_threads) {
              // NB: we will write data as threads process the data.  This is to
              // facilitate SAM/BAM writing, which may be slow, especially for
              // BAM.
              int32_t tid = (i % driver->opt->num_threads);
              while((*thread_data[tid].buffer_idx) <= i) {
                  usleep(1000*1000); // sleep
              }
          }
#endif
          // write
          for(j=0;j<bams[i]->n;j++) { // for each end
              for(k=0;k<bams[i]->bams[j]->n;k++) { // for each hit
                  bam1_t *b = NULL;
                  b = bams[i]->bams[j]->bams[k]; // that's a lot of BAMs
                  if(samwrite(io_out->fp, b) <= 0) {
                      tmap_error("Error writing the SAM file", Exit, WriteFileError);
                  }
              }
          }
          tmap_map_bams_destroy(bams[i]);
          bams[i] = NULL;
      }

#ifdef HAVE_LIBPTHREAD
      // join threads
      if(1 <driver->opt->num_threads) {
          // join threads
          for(i=0;i<driver->opt->num_threads;i++) {
              if(0 != pthread_join(threads[i], NULL)) {
                  tmap_error("error joining threads", Exit, ThreadError);
              }
              // add the stats
              tmap_map_stats_add(stat, stats[i]);
              // free the buffer index
              free(thread_data[i].buffer_idx);
          }
          free(threads); threads = NULL;
          free(thread_data); thread_data = NULL;
          free(attr); attr = NULL;
      }
#endif
      // TODO: should we flush when writing SAM and processing one read at a time?

      // print statistics
      n_reads_processed += seqs_buffer_length;
      if(-1 != driver->opt->reads_queue_size) {
          tmap_progress_print2("processed %d reads", n_reads_processed);
          tmap_progress_print2("stats [%.2lf,%.2lf,%.2lf,%.2lf,%.2lf]",
                               stat->num_with_mapping * 100.0 / (double)stat->num_reads,
                               stat->num_after_seeding/(double)stat->num_with_mapping,
                               stat->num_after_scoring/(double)stat->num_with_mapping,
                               stat->num_after_rmdup/(double)stat->num_with_mapping,
                               stat->num_after_filter/(double)stat->num_with_mapping);
      }
  }
  if(-1 == driver->opt->reads_queue_size) {
      tmap_progress_print2("processed %d reads", n_reads_processed);
      tmap_progress_print2("stats [%.2lf,%.2lf,%.2lf,%.2lf,%.2lf]",
                           stat->num_with_mapping * 100.0 / (double)stat->num_reads,
                           stat->num_after_seeding/(double)stat->num_with_mapping,
                           stat->num_after_scoring/(double)stat->num_with_mapping,
                           stat->num_after_rmdup/(double)stat->num_with_mapping,
                           stat->num_after_filter/(double)stat->num_with_mapping);
  }

  // cleanup the algorithm persistent data
  tmap_map_driver_do_cleanup(driver);

  // close the input/output
  tmap_sam_io_destroy(io_out);

  // free memory
  tmap_index_destroy(index);
  tmap_seqs_io_destroy(io_in);
  for(i=0;i<reads_queue_size;i++) {
      tmap_seqs_destroy(seqs_buffer[i]);
  }
  free(seqs_buffer);
  free(records);
  free(bams);
  tmap_map_stats_destroy(stat);
#ifdef HAVE_LIBPTHREAD
  for(i=0;i<driver->opt->num_threads;i++) {
      tmap_map_stats_destroy(stats[i]);
      tmap_rand_destroy(rand[i]);
  }
  free(stats);
  free(rand);
#else
  tmap_rand_destroy(rand);
#endif
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
  tmap_rand_destroy(rand_core);
#endif
}

/* MAIN API */

tmap_map_driver_algorithm_t*
tmap_map_driver_algorithm_init(tmap_map_driver_func_init func_init,
                    tmap_map_driver_func_thread_init func_thread_init,
                    tmap_map_driver_func_thread_map func_thread_map,
                    tmap_map_driver_func_thread_cleanup func_thread_cleanup,
                    tmap_map_driver_func_cleanup func_cleanup,
                    tmap_map_opt_t *opt)
{
  tmap_map_driver_algorithm_t *algorithm = NULL;

  // the only necessary function is func_thread_map
  if(func_thread_map == NULL ) {
      tmap_error("func_thread_map cannot be null", Exit, OutOfRange);
  }

  algorithm = tmap_calloc(1, sizeof(tmap_map_driver_algorithm_t), "algorithm");
  algorithm->func_init = func_init;
  algorithm->func_thread_init = func_thread_init;
  algorithm->func_thread_map = func_thread_map;
  algorithm->func_thread_cleanup = func_thread_cleanup;
  algorithm->func_cleanup = func_cleanup;
  algorithm->opt = opt;
  algorithm->data = NULL;
  algorithm->thread_data = tmap_calloc(opt->num_threads, sizeof(void*), "algorithm->thread_data");
  return algorithm;
}

void
tmap_map_driver_algorithm_destroy(tmap_map_driver_algorithm_t *algorithm)
{
  free(algorithm->thread_data);
  free(algorithm);
}

tmap_map_driver_stage_t*
tmap_map_driver_stage_init(int32_t stage)
{
  tmap_map_driver_stage_t *s = NULL;
  s = tmap_calloc(1, sizeof(tmap_map_driver_stage_t), "stage");
  s->stage = stage;
  s->opt = tmap_map_opt_init(TMAP_MAP_ALGO_STAGE);
  s->opt->algo_stage = stage;
  return s;
}

void
tmap_map_driver_stage_add(tmap_map_driver_stage_t *s,
                    tmap_map_driver_func_init func_init,
                    tmap_map_driver_func_thread_init func_thread_init,
                    tmap_map_driver_func_thread_map func_thread_map,
                    tmap_map_driver_func_thread_cleanup func_thread_cleanup,
                    tmap_map_driver_func_cleanup func_cleanup,
                    tmap_map_opt_t *opt)
{
  // check against stage options
  tmap_map_opt_check_stage(s->opt, opt);
  s->num_algorithms++;
  s->algorithms = tmap_realloc(s->algorithms, sizeof(tmap_map_driver_algorithm_t*) * s->num_algorithms, "s->algorithms");
  s->algorithms[s->num_algorithms-1] = tmap_map_driver_algorithm_init(func_init, func_thread_init, func_thread_map,
                                                                      func_thread_cleanup, func_cleanup, opt);
}

void
tmap_map_driver_stage_destroy(tmap_map_driver_stage_t *stage)
{
  int32_t i;
  for(i=0;i<stage->num_algorithms;i++) {
      tmap_map_driver_algorithm_destroy(stage->algorithms[i]);
  }
  tmap_map_opt_destroy(stage->opt);
  free(stage->algorithms);
  free(stage);
}

tmap_map_driver_t*
tmap_map_driver_init(int32_t algo_id, tmap_map_driver_func_mapq func_mapq)
{
  tmap_map_driver_t *driver = NULL;
  driver = tmap_calloc(1, sizeof(tmap_map_driver_t), "driver");
  driver->func_mapq = func_mapq;
  driver->opt = tmap_map_opt_init(algo_id);
  return driver;
}

void
tmap_map_driver_add(tmap_map_driver_t *driver,
                    tmap_map_driver_func_init func_init,
                    tmap_map_driver_func_thread_init func_thread_init,
                    tmap_map_driver_func_thread_map func_thread_map,
                    tmap_map_driver_func_thread_cleanup func_thread_cleanup,
                    tmap_map_driver_func_cleanup func_cleanup,
                    tmap_map_opt_t *opt)
{
  // make more stages
  if(driver->num_stages < opt->algo_stage) {
      driver->stages = tmap_realloc(driver->stages, sizeof(tmap_map_driver_stage_t*) * opt->algo_stage, "driver->stages");
      while(driver->num_stages < opt->algo_stage) {
          driver->num_stages++;
          driver->stages[driver->num_stages-1] = tmap_map_driver_stage_init(driver->num_stages);
          // copy global options into this stage
          tmap_map_opt_copy_global(driver->stages[driver->num_stages-1]->opt, driver->opt);
          // copy stage options into this stage
          tmap_map_opt_copy_stage(driver->stages[driver->num_stages-1]->opt, opt);
      }
  }

  // check options
  tmap_map_opt_check(opt);

  // check against global options
  tmap_map_opt_check_global(driver->opt, opt);

  // add to the stage
  tmap_map_driver_stage_add(driver->stages[opt->algo_stage-1],
                    func_init,
                    func_thread_init,
                    func_thread_map,
                    func_thread_cleanup,
                    func_cleanup,
                    opt);
}

void
tmap_map_driver_run(tmap_map_driver_t *driver)
{
  tmap_map_driver_core(driver);
}

void
tmap_map_driver_destroy(tmap_map_driver_t *driver)
{
  int32_t i;
  for(i=0;i<driver->num_stages;i++) {
      tmap_map_driver_stage_destroy(driver->stages[i]);
  }
  tmap_map_opt_destroy(driver->opt);
  free(driver->stages);
  free(driver);
}
