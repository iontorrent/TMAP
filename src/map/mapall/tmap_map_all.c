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
#include <limits.h>
#include "../../util/tmap_error.h"
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_definitions.h"
#include "../../util/tmap_progress.h"
#include "../../util/tmap_sam_print.h"
#include "../../util/tmap_sort.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../../index/tmap_bwt_gen.h"
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_bwt_match.h"
#include "../../index/tmap_bwt_match_hash.h"
#include "../../index/tmap_sa.h"
#include "../../index/tmap_index.h"
#include "../../io/tmap_seq_io.h"
#include "../../server/tmap_shm.h"
#include "../../sw/tmap_sw.h"
#include "../util/tmap_map_stats.h"
#include "../util/tmap_map_util.h"
#include "../map1/tmap_map1.h"
#include "../map1/tmap_map1_aux.h"
#include "../map2/tmap_map2.h"
#include "../map2/tmap_map2_aux.h"
#include "../map3/tmap_map3.h"
#include "../map3/tmap_map3_aux.h"
#include "../mapvsw/tmap_map_vsw.h"
#include "../tmap_map_driver.h"
#include "tmap_map_all.h"

#define tmap_map_all_stages_used_lt(a, b) (a < b) 
TMAP_SORT_INIT(tmap_map_all_stages_used, int32_t, tmap_map_all_stages_used_lt)

static void
tmap_map_all_add_algorithm(tmap_map_driver_t *driver, tmap_map_opt_t *opt)
{
  switch(opt->algo_id) {
    case TMAP_MAP_ALGO_MAP1:
      // add this algorithm
      tmap_map_driver_add(driver,
                          tmap_map1_init,
                          tmap_map1_thread_init,
                          tmap_map1_thread_map,
                          tmap_map1_thread_cleanup,
                          NULL,
                          opt);
      break;
    case TMAP_MAP_ALGO_MAP2:
      // add this algorithm
      tmap_map_driver_add(driver,
                          tmap_map2_init,
                          tmap_map2_thread_init,
                          tmap_map2_thread_map,
                          tmap_map2_thread_cleanup,
                          NULL,
                          opt);
      break;
    case TMAP_MAP_ALGO_MAP3:
      // add this algorithm
      tmap_map_driver_add(driver,
                          tmap_map3_init,
                          tmap_map3_thread_init,
                          tmap_map3_thread_map,
                          tmap_map3_thread_cleanup,
                          NULL,
                          opt);
      break;
    case TMAP_MAP_ALGO_MAPVSW:
      // add this algorithm
      tmap_map_driver_add(driver,
                          tmap_map_vsw_init,
                          tmap_map_vsw_thread_init,
                          tmap_map_vsw_thread_map,
                          tmap_map_vsw_thread_cleanup,
                          NULL,
                          opt);
      break;
    case TMAP_MAP_ALGO_STAGE:
      // ignore
      break;
    default:
      tmap_error("Unknown algorithm", Exit, OutOfRange);
  }
}

static void
tmap_map_all_core(tmap_map_driver_t *driver)
{
  int32_t i;

  // add the algorithms
  for(i=0;i<driver->opt->num_sub_opts;i++) {
      tmap_map_all_add_algorithm(driver, driver->opt->sub_opts[i]);
  }

  // DEBUG
  /*
  int32_t j;
  for(i=0;i<driver->num_stages;i++) {
      for(j=0;j<driver->stages[i]->num_algorithms;j++) {
          fprintf(stderr, "Algorithm: %s Stage: %d\n", 
                  tmap_algo_id_to_name(driver->stages[i]->algorithms[j]->opt->algo_id),
                  driver->stages[i]->algorithms[j]->opt->algo_stage);
      }
  }
  */

  // run the driver
  tmap_map_driver_run(driver);
}

static int32_t
tmap_map_all_opt_get_next_stage_j(int argc, char *argv[], int32_t j)
{
  int32_t cur_stage;
  while(j < argc) {
      if(6 <= strlen(argv[j]) 
         && 0 == strncmp("stage", argv[j], 5)) {
          cur_stage = atoi(argv[j] + 5);
          // get the stage name
          if(cur_stage <= 0 || INT_MAX == cur_stage) {
              tmap_error("Could not identify the stage", Exit, CommandLineArgument);
          }
          if(strlen(argv[j]) == 6 + (int)log10(cur_stage)) { 
              break;
          }
      }
      j++;
  }
  return j;
}

int32_t
tmap_map_all_opt_parse(int argc, char *argv[], tmap_map_opt_t *opt)
{
  int32_t i, j, k, l;
  int32_t cur_id, cur_stage;
  int32_t *stages_used = NULL, stages_used_num = 0;
  tmap_map_opt_t *stage_opt = NULL, *algo_opt = NULL;

  // DEBUG
  /*
  for(i=0;i<argc;i++) {
      fprintf(stderr, "argv[%d]=%s\n", i, argv[i]);
  }
  */

  // get the range for the global options
  i = 1;
  j = tmap_map_all_opt_get_next_stage_j(argc, argv, 1);
  if(j == argc) {
      tmap_error("No stages were specified", Warn, CommandLineArgument);
  }

  // parse the global options
  optind = 1; 
  if(0 == tmap_map_opt_parse(j, argv, opt)) {
      return 0;
  }

  // go through the stages
  i = j;
  while(i < argc) {
      // find the next stage, if it exists
      j = tmap_map_all_opt_get_next_stage_j(argc, argv, i+1);
          
      // get the stage name
      cur_stage = atoi(argv[i] + 5);
      if(cur_stage <= 0 || INT_MAX == cur_stage) {
          tmap_error("Could not identify the stage", Exit, CommandLineArgument);
      }

      // check that the stage was not previously specified
      for(k=0;k<stages_used_num;k++) {
          if(cur_stage == stages_used[k]) {
              tmap_error("Cannot specify the same stage twice", Exit, CommandLineArgument);
          }
      }

      // add to the stages used
      stages_used_num++;
      stages_used = tmap_realloc(stages_used, sizeof(uint32_t) * stages_used_num, "stages_used");
      stages_used[stages_used_num-1] = cur_stage;
      
      // get the stage options
      k = i+1; // ignore stage name
      while(k < j) {
          cur_id = tmap_algo_name_to_id(argv[k]);  // check for an algorithm id
          if(0 < cur_id) break;
          k++;
      }

      // parse the stage options
      // from i to k
      stage_opt = tmap_map_opt_add_sub_opt(opt, TMAP_MAP_ALGO_STAGE);
      stage_opt->algo_stage = cur_stage;
      optind = 1; 
      if(1 != tmap_map_opt_parse(k - i, argv + i, stage_opt)) {
          return 0;
      }

      // NB: do this after parsing the stage options, since '-h' may be
      // specified
      if(k == j) {
          tmap_error("A stage was specified with no algorithms", Exit, CommandLineArgument);
      }
      
      // parse the algorithms in this stage
      while(k < j) {
          // get the algorithm id
          cur_id = tmap_algo_name_to_id(argv[k]); 
          if(cur_id <= 0) tmap_error("bug encountered", Exit, OutOfRange); // should not happen
          algo_opt = tmap_map_opt_add_sub_opt(opt, cur_id);
          algo_opt->algo_stage = cur_stage;

          // copy global options to the algorithm options
          tmap_map_opt_copy_stage(algo_opt, opt);
          // copy stage options to the algorithm options
          tmap_map_opt_copy_stage(algo_opt, stage_opt);
          
          // check that the algorithm has not been specified
          for(l=0;l<opt->num_sub_opts;l++) {
              if(opt->sub_opts[l]->algo_id == algo_opt->algo_id 
                 && opt->sub_opts[l]->algo_stage == opt->algo_stage) {
                  tmap_error("algorithm specified twice for the same stage", Exit, CommandLineArgument);
              }
          }

          // get the range of options
          l = k+1;
          while(l < j && tmap_algo_name_to_id(argv[l]) <= 0) {
              l++;
          }

          // parse the algorithm options
          optind = 1; 
          if(1 != tmap_map_opt_parse(l - k, argv + k, algo_opt)) {
              return 0;
          }

          // for the next loop
          k = l;

          // nullify
          algo_opt = NULL;
      }

      // nullify
      stage_opt = NULL;

      // update the outer loop
      i = j;
  }

  // sort stages used
  tmap_sort_introsort(tmap_map_all_stages_used, stages_used_num, stages_used);

  // check that stages are continuous
  for(i=0;i<stages_used_num;i++) {
      if(i+1 != stages_used[i]) {
          tmap_error("stage was missing", Exit, CommandLineArgument);
      }
  }

  // do this after parsing
  opt->argc = argc; opt->argv = argv;

  // free
  free(stages_used);

  optind = argc;
  return 1;
}

int 
tmap_map_all_main(int argc, char *argv[])
{
  tmap_map_driver_t *driver = NULL;

  // init opt
  driver = tmap_map_driver_init(TMAP_MAP_ALGO_MAPALL, tmap_map_util_mapq);

  // get options
  if(1 != tmap_map_all_opt_parse(argc, argv, driver->opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return tmap_map_opt_usage(driver->opt);
  }
  else { 
      // check command line arguments
      tmap_map_opt_check(driver->opt);
  }

  // run map_all
  tmap_map_all_core(driver);

  // destroy 
  tmap_map_driver_destroy(driver);

  tmap_progress_print2("terminating successfully");

  return 0;
}
