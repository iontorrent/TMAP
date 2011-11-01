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
#include "../../util/tmap_error.h"
#include "../../util/tmap_alloc.h"
#include "../../util/tmap_definitions.h"
#include "../../util/tmap_progress.h"
#include "../../util/tmap_sam_print.h"
#include "../../seq/tmap_seq.h"
#include "../../index/tmap_refseq.h"
#include "../../index/tmap_bwt_gen.h"
#include "../../index/tmap_bwt.h"
#include "../../index/tmap_bwt_match.h"
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
                          tmap_map_util_mapq,
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
                          tmap_map_util_mapq,
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
                          tmap_map_util_mapq,
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
                          tmap_map_util_mapq,
                          tmap_map_vsw_thread_cleanup,
                          NULL,
                          opt);
      break;
    default:
      tmap_error("Unknown algorithm", Exit, OutOfRange);
  }
}

static void
tmap_map_all_core(tmap_map_opt_t *opt)
{
  int32_t i;
  tmap_map_driver_t *driver = NULL;

  driver = tmap_map_driver_init();

  // add the algorithms
  for(i=0;i<opt->num_sub_opts;i++) {
      tmap_map_all_add_algorithm(driver, opt->sub_opts[i]);
  }

  // run the driver
  tmap_map_driver_run(driver);

  tmap_map_driver_destroy(driver);
}

int32_t
tmap_map_all_opt_parse(int argc, char *argv[], tmap_map_opt_t *opt)
{
  int32_t i, j, start, opt_id, opt_next_id, opt_stage, opt_next_stage;
  char *name = NULL;
  tmap_map_opt_t *opt_next = NULL;

  // parse common options as well as map1/map2/map3/mapvsw commands
  start = 0;
  i = 1;
  opt_id = opt_next_id = TMAP_MAP_ALGO_NONE;
  opt_stage = opt_next_stage = 0;
  opt->num_stages = 0;
  while(i<=argc) {
      if(i == argc) { 
          // do nothing
      }
      else {
          // get the algorithm type and stage
          opt_next_stage = 1;
          name = tmap_strdup(argv[i]); // copy the command line option
          while(opt_next_stage <= 2) { // while it could be the first or second stage
              j = tmap_algo_name_to_id(name); // get the algorithm id
              if(0 < j) { // found!
                  opt_next_id = j;
                  if(opt_next_stage-1 == opt->num_stages) opt->num_stages++;
                  break;
              }
              if(1 == opt_next_stage) {
                  // convert to lower case
                  for(j=0;j<strlen(name);j++) {
                      name[j] = tolower(name[j]);
                  }
              }
              opt_next_stage++;
          }
          if(2 < opt_next_stage) opt_next_stage = 0;
          free(name);
      }
      
      /*
      fprintf(stderr, "i=%d start=%d argc=%d opt_id=%d name=%s opt_next_id=%d name=%s argv[%d]=%s argc=%d\n",
              i, start, argc, 
              opt_id, tmap_algo_id_to_name(opt_id),
              opt_next_id, tmap_algo_id_to_name(opt_next_id), 
              i, argv[i], argc);
              */
      
      if(opt_id != opt_next_id // new type
         || opt_stage != opt_next_stage // new stage
         || i == argc) { // end of command line arguments
          optind=1; // needed for getopt_long

          /*
          int j;
          for(j=0;j<i-start;j++) {
              fprintf(stderr, "j=%d arg=%s\n", j, argv[j+start]);
          }
          */

          if(TMAP_MAP_ALGO_NONE == opt_id) {
              // parse common options
              if(0 == tmap_map_opt_parse(i-start, argv+start, opt)) {
                  return 0;
              }
          }
          else {
              // get a sub-opt
              opt_next = tmap_map_opt_add_sub_opt(opt, opt_next_id);
              // set id and stage
              opt_next->algo_id = opt_next_id;
              opt_next->algo_stage = opt_next_stage;
              // parse common options
              if(0 < i - start) {
                  if(0 == tmap_map_opt_parse(i-start, argv+start, opt_next)) {
                      return 0;
                  }
              }
          }

          // update next
          opt_id = opt_next_id;
          opt_stage = opt_next_stage;
          start = i;
      }
      i++;
  }
  if(argc < i) {
      i = argc;
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
