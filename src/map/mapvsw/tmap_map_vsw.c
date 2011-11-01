/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <math.h>
#include <config.h>
#ifdef HAVE_LIBPTHREAD
#include <pthread.h>
#endif
#include <unistd.h>
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
#include "../../index/tmap_sa.h"
#include "../../index/tmap_index.h"
#include "../../io/tmap_seq_io.h"
#include "../../server/tmap_shm.h"
#include "../../sw/tmap_sw.h"
#include "../util/tmap_map_stats.h"
#include "../util/tmap_map_util.h"
#include "../tmap_map_driver.h"
#include "tmap_map_vsw.h"

static int32_t
tmap_map_vsw_mapq(tmap_map_sams_t *sams, int32_t seq_len, tmap_map_opt_t *opt)
{
  tmap_map_util_mapq(sams, seq_len, opt);
  return 0;
}

int32_t
tmap_map_vsw_init(void **data, tmap_refseq_t *refseq, tmap_map_opt_t *opt)
{
  return 0;
}

int32_t 
tmap_map_vsw_thread_init(void **data, tmap_map_opt_t *opt)
{
  (*data) = NULL;
  return 0;
}

// reverse and reverse compliment
tmap_map_sams_t*
tmap_map_vsw_thread_map_core(void **data, tmap_seq_t *seqs[2], int32_t seq_len,
                             tmap_index_t *index, tmap_map_opt_t *opt)
{
  int32_t i;
  tmap_map_sams_t *sams = NULL;

  if((0 < opt->min_seq_len && seq_len < opt->min_seq_len)
     || (0 < opt->max_seq_len && opt->max_seq_len < seq_len)) {
      // go to the next loop
      return tmap_map_sams_init(NULL);
  }

  sams = tmap_map_sams_init(NULL);
  tmap_map_sams_realloc(sams, index->refseq->num_annos<<1); // one for each contig and strand

  for(i=0;i<index->refseq->num_annos<<1;i++) {
      tmap_map_sam_t *s;

      // save
      s = &sams->sams[i];

      // save the hit
      s->algo_id = TMAP_MAP_ALGO_MAPVSW;
      s->algo_stage = 0;
      s->strand = i & 1;
      s->seqid = i >> 1;
      s->pos = 0;
      s->target_len = index->refseq->annos[s->seqid].len;
      s->score_subo = INT32_MIN;
      
      // mapvswaux data
      tmap_map_sam_malloc_aux(s, TMAP_MAP_ALGO_MAPVSW);
  }

  return sams;
}

tmap_map_sams_t*
tmap_map_vsw_thread_map(void **data, tmap_seq_t **seqs, tmap_index_t *index, tmap_map_stats_t *stat, tmap_rand_t *rand, tmap_map_opt_t *opt)
{
  int32_t seq_len = 0;;
  tmap_seq_t *seqs_tmp[2]={NULL, NULL};
  tmap_map_sams_t *sams = NULL;

  // sequence length
  seq_len = tmap_seq_get_bases(seqs[0])->l;

  // sequence length not in range
  if((0 < opt->min_seq_len && seq_len < opt->min_seq_len)
     || (0 < opt->max_seq_len && opt->max_seq_len < seq_len)) {
      return tmap_map_sams_init(NULL);
  }

  // clone the sequence 
  seqs_tmp[0] = seqs[0];
  seqs_tmp[1] = seqs[1];

  // core algorithm
  sams = tmap_map_vsw_thread_map_core(data, seqs_tmp, seq_len, index, opt);

  return sams;
}

int32_t
tmap_map_vsw_thread_cleanup(void **data, tmap_map_opt_t *opt)
{
  return 0;
}

static void 
tmap_map_vsw_core(tmap_map_opt_t *opt)
{
  tmap_map_driver_t *driver = NULL;

  driver = tmap_map_driver_init();

  // add this algorithm
  tmap_map_driver_add(driver,
                      tmap_map_vsw_init, 
                      tmap_map_vsw_thread_init, 
                      tmap_map_vsw_thread_map, 
                      tmap_map_vsw_mapq,
                      tmap_map_vsw_thread_cleanup,
                      NULL,
                      opt);
  

  // run the driver
  tmap_map_driver_run(driver);

  tmap_map_driver_destroy(driver);
}

int 
tmap_map_vsw_main(int argc, char *argv[])
{
  tmap_map_opt_t *opt = NULL;

  // init opt
  opt = tmap_map_opt_init(TMAP_MAP_ALGO_MAPVSW);

  // get options
  if(1 != tmap_map_opt_parse(argc, argv, opt) // options parsed successfully
     || argc != optind  // all options should be used
     || 1 == argc) { // some options should be specified
      return tmap_map_opt_usage(opt);
  }
  else { 
      // check command line arguments
      tmap_map_opt_check(opt);
  }

  // run map_vsw
  tmap_map_vsw_core(opt);

  // destroy opt
  tmap_map_opt_destroy(opt);

  tmap_progress_print2("terminating successfully");

  return 0;
}
