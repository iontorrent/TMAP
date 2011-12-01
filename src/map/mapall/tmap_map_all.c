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

  /*
  for(i=0;i<driver->num_algorithms;i++) {
      fprintf(stderr, "Algorithm: %s Stage: %d\n", 
              tmap_algo_id_to_name(driver->algorithms[i]->opt->algo_id),
              driver->algorithms[i]->opt->algo_stage);
  }
  */

  // run the driver
  tmap_map_driver_run(driver);
}

void
tmap_map_all_default_opts(tmap_map_opt_t* opt) {
    /*
     * int32_t num_stages;  !< the number of stages 
    int32_t mapall_score_thr;  !< the stage one scoring threshold (match-score-scaled) (--staged-score-thres) 
    int32_t mapall_mapq_thr;  !< the stage one mapping quality threshold (--staged-mapq-thres) 
    int32_t mapall_keep_all;  !< keep mappings that do not pass the first stage threshold for the next stage (--staged-keep-all) 
    double  mapall_seed_freqc; !< the minimum frequency a seed must occur in order to be considered for mapping (--seed-freq-cutoff) 
     */
    
    opt->mapall_score_thr = 8;
    opt->mapall_mapq_thr = 23;
    opt->mapall_keep_all = 1;
    opt->mapall_seed_freqc = 0.0;
}

void
tmap_map_all_copy_stage_opts(tmap_map_opt_t* src, tmap_map_opt_t* dest) {
    src->algo_stage = dest->algo_stage;
    src->algo_id = dest->algo_id;
    dest->mapall_score_thr = src->mapall_score_thr;
    dest->mapall_mapq_thr = src->mapall_mapq_thr;
    dest->mapall_keep_all = src->mapall_keep_all;
    dest->mapall_seed_freqc = src->mapall_seed_freqc;
    
}

int32_t
tmap_map_all_opt_parse(int argc, char *argv[], tmap_map_opt_t *opt)
{
  int32_t i, j, z, start, opt_id, opt_next_id, opt_stage, opt_next_stage, cur_id, cur_stage;
  char *name = NULL;
  tmap_map_opt_t *opt_cur= NULL;
  //refactor options
  int32_t n_stages, n, stage_starts_index, stage_stops_index;
  int32_t* stage_starts;
  int32_t* stage_stops;
  int32_t mapall_id = tmap_algo_name_to_id("mapall");
  tmap_map_opt_t *mapall_ops; 
  // parse common options as well as map1/map2/map3/mapvsw commands
  start = 0;
  i = 1;
  opt_id = opt_next_id = TMAP_MAP_ALGO_NONE;
  opt_stage = opt_next_stage = 0;
  opt->num_stages = 0;
  
  n = n_stages = stage_starts_index = stage_stops_index = 0;
  while(n < argc) {
      //if(0 == strcmp(name, algo_id_to_name[i])) return id;

      if (0 == strcmp("stage", argv[n])) {
          n_stages++;
      }
      n++;
  }
  printf("n: %d n_stages: %d\n", n, n_stages);
  if (n_stages > 0) {
      //setup arrays
      stage_starts = tmap_malloc(n_stages * sizeof(int32_t), "stage_starts" );
      stage_stops = tmap_malloc(n_stages * sizeof(int32_t), "stage_stops" );
      
      /*get start and stop positions in argv array 
        of stage locations*/
      for (n = 0; n < argc; n++) {
        if (0 == strcmp("stage", argv[n])) {
          stage_starts[ stage_starts_index++ ] = n;

          if (stage_starts_index > 1) {
              printf("stop_index: %d start_index: %d\n", stage_stops_index, stage_starts_index);
              stage_stops[ stage_stops_index++ ] = stage_starts[ stage_starts_index - 1] - 1;
          }

        }       
      }
      stage_stops[ stage_stops_index ] = n - 1;
      
      
  } else { // what the hell do i do here?
      
  }
  
  for (n=0; n < n_stages; n++) {
      printf("stage: %d  start: %d  stop: %d argc: %d\n\t", n, stage_starts[n], stage_stops[n], argc);
      for(i=stage_starts[n]; i <= stage_stops[n]; i++){
          printf("%s ", argv[i]);
      }
      printf("\n");
  }
  cur_stage = 0;
  //parse common options
  printf("common opts:\n\t");
  for (n=1; n < stage_starts[0]; n++){
      printf("%s ", argv[n]);
  }
  printf("\n");
 // printf("1 - stage_starts[0]: %d\n", stage_starts[0] - 1);
  if(0 == tmap_map_opt_parse(stage_starts[0], argv, opt)) {
      //no common options in between mapall and stage 1, or 
      //possibly malformed options all together
        return 0;
  }
  

  
  for(z=0; z < n_stages; z++) { //z is the stage
      printf("z: %d ", z);
      mapall_ops = tmap_malloc(sizeof(tmap_map_opt_t), "mapall_ops");
      tmap_map_all_default_opts(mapall_ops);
      cur_stage = atoi(argv[stage_starts[z]+1]);
      printf("n_stages: %d argv: %s ", n_stages, argv[stage_starts[z]+1]);
      printf("cur_stage: %d stage_starts[%d]=%s  stage_stops[%d]=%s\n", 
              cur_stage, z, argv[stage_starts[z]], z, argv[stage_stops[z]]);
      for(i = stage_starts[z] + 2; i <= stage_stops[z]; i++) {

          // get the algorithm type and stage
          name = tmap_strdup(argv[i]); // copy the command line option
          
          cur_id = tmap_algo_name_to_id(name); // get the algorithm id
          printf("cur_id: %d argv[%d]: %s\n", cur_id, i, argv[i] );

          //if it's a mapping algo
          if(0 < cur_id) { // found!
              start = i;
              //printf("cur_id: %d argv[start:%d]: %s\n", cur_id, start, argv[start]);
              opt_id = cur_id;
              i++;//move i to next element
              while ( i <= stage_stops[z] && 0 > tmap_algo_name_to_id( argv[i] ) ) {
                  printf("if: argv[i:%d] %s\n", i, argv[i]);
                  i++;                  
              }
              
              
          }//we're not on an algorithm
          else {
              start = i;
              //printf("else argv[start:%d]: %s\n", start, argv[start]);
              //opt_id = mapall_id;
              
              while( i <= stage_stops[z] && 0 > tmap_algo_name_to_id( argv[i] )  ) {
                  printf("else: argv[i:%d] %s\n", i, argv[i]);
                  i++;
              }
              //rewind to the map algo name
              printf("i: %d  stage_stops[z]: %d\n", i, stage_stops[z]);
              if (i !=stage_stops[z] + 1){
                  i--;
                  printf("rewind!\n");
              }
              continue; //next iteration of for loop
          }
          
          free(name);
      
      /*
      fprintf(stderr, "ITER i=%d start=%d argc=%d opt_id=%d name=%s opt_stage=%d opt_next_id=%d name=%s opt_next_stage=%d argv[%d]=%s argc=%d\n",
              i, start, argc, 
              opt_id, tmap_algo_id_to_name(opt_id), opt_stage,
              opt_next_id, tmap_algo_id_to_name(opt_next_id), opt_next_stage, 
              i, (i < argc) ? argv[i] : NULL, argc);
      */
      
      /*if(opt_id != opt_next_id // new type
         || i == argc) { // end of command line arguments*/
          optind=1; // needed for getopt_long
          /*
          fprintf(stderr, "Algorithm: %s start=%d i=%d\n", tmap_algo_id_to_name(opt_id), start, i);
          int j;
          for(j=0;j<i-start;j++) {
              fprintf(stderr, "j=%d arg=%s\n", j, argv[j+start]);
          }
          */
          if (opt_id != mapall_id) {
          printf("in lower if\n");
          if(opt->num_stages < cur_stage) { 
              opt->num_stages = cur_stage;
              
          }

          if(0 < i - start) {
              if (opt_id == mapall_id) {
                  //optind=0;
                  //printf("argv+start=%s \n", *(argv+start) );

                  if(0 == tmap_map_opt_parse((i+1)-start, argv+start-1, mapall_ops)) {
                    return 0;
                  }
              }
              else{

                // printf("algo_id: %d argv[i:%d]: %s  argv[start:%d]: %s  i - start: %d\n", opt_cur->algo_id,i,argv[i], start,argv[start], i - start);
                  printf("argv[lawl i: %d] %s\n", i, argv[i]);
                  opt_cur= tmap_map_opt_add_sub_opt(opt, opt_id);

                  // set stage
                  opt_cur->algo_stage = cur_stage;
                  //set algo id
                  opt_cur->algo_id = opt_id;
                  // parse common options
                  if( 0 == tmap_map_opt_parse(i-start, argv+start, opt_cur) ) {
                    return 0;
                  } else {
                      tmap_map_all_copy_stage_opts(mapall_ops, opt_cur);
                  }
              }
          }


          i--;
      }

    }//end for over stage start/stop range

  }//end for over n_stages
  if(argc < i) {
      i = argc;
  }
  optind = i;
  
  // do this after parsing
  opt->argc = argc; opt->argv = argv;
  //for i = num_sub_opts; 
   //sub_opts[i]->algo_id;
  
  for(i=0; i<opt->num_sub_opts; i++){
      printf("sub opt algo_id: %d stage: %d\n", opt->sub_opts[i]->algo_id, opt->sub_opts[i]->algo_stage);
  }
  //free(mapall_ops);
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
