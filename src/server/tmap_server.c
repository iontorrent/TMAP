#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <unistd.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>

#include "../util/tmap_error.h"
#include "../util/tmap_alloc.h"
#include "../util/tmap_progress.h"
#include "../index/tmap_refseq.h"
#include "../index/tmap_bwt.h"
#include "../index/tmap_sa.h"
#include "tmap_shm.h"
#include "tmap_server.h"

static tmap_shm_t *tmap_server_shm_ptr = NULL;

void
tmap_server_sigint(int signal)
{
  if(NULL != tmap_server_shm_ptr) {
      tmap_error("will try to destroy the shared memory", Warn, SigInt); // warn
      // try to destroy the shared memory
      tmap_shm_set_dead(tmap_server_shm_ptr);
      tmap_shm_destroy(tmap_server_shm_ptr, 1);
      tmap_server_shm_ptr = NULL;
      tmap_error("shared memory destroyed", Exit, SigInt); //exit
  }
  else {
      tmap_error(NULL, Exit, SigInt); //exit
  }
}

void
tmap_server_set_sigint(tmap_shm_t *shm)
{
  tmap_server_shm_ptr = shm;
  signal(SIGINT, tmap_server_sigint); 
}

void
tmap_server_start(char *fn_fasta, key_t key, uint32_t listing)
{
  int32_t i;
  tmap_shm_t *shm = NULL;
  tmap_refseq_t *refseq = NULL;
  tmap_bwt_t *bwt = NULL;
  tmap_sa_t *sa = NULL;
  uint8_t *buf = NULL;
  size_t n_bytes = 0, cur_bytes = 0;
  uint32_t cur_listing = 0;

  tmap_progress_print("starting server");

  // get data size
  n_bytes = 0;
  if(listing & TMAP_SHM_LISTING_REFSEQ) {
      n_bytes += tmap_refseq_shm_read_num_bytes(fn_fasta, 0);
  }
  if(listing & TMAP_SHM_LISTING_REV_REFSEQ) {
      n_bytes += tmap_refseq_shm_read_num_bytes(fn_fasta, 1);
  }
  if(listing & TMAP_SHM_LISTING_BWT) {
      n_bytes += tmap_bwt_shm_read_num_bytes(fn_fasta, 0);
  }
  if(listing & TMAP_SHM_LISTING_REV_BWT) {
      n_bytes += tmap_bwt_shm_read_num_bytes(fn_fasta, 1);
  }
  if(listing & TMAP_SHM_LISTING_SA) {
      n_bytes += tmap_sa_shm_read_num_bytes(fn_fasta, 0);
  }
  if(listing & TMAP_SHM_LISTING_REV_SA) {
      n_bytes += tmap_sa_shm_read_num_bytes(fn_fasta, 1);
  }

  // get shared memory
  tmap_progress_print("retrieving shared memory [%llu bytes]", (long long unsigned int)n_bytes);
  shm = tmap_shm_init(key, n_bytes, 1);
  tmap_progress_print2("shared memory retrieved");

  // catch a ctrl-c signal from now on
  tmap_server_set_sigint(shm);

  // pack the shared memory
  tmap_progress_print("packing shared memory");
  buf = (uint8_t*)shm->buf;

  // pack the reference sequence
  for(i=0;i<2;i++) { // forward/reverse
      cur_listing = (0 == i) ? TMAP_SHM_LISTING_REFSEQ : TMAP_SHM_LISTING_REV_REFSEQ;
      if(listing & cur_listing) {
          tmap_progress_print("packing %s reference", (0 == i) ? "forward" : "reverse");
          refseq = tmap_refseq_read(fn_fasta, i);
          cur_bytes = tmap_refseq_shm_num_bytes(refseq);
          tmap_refseq_shm_pack(refseq, buf);
          tmap_shm_add_listing(shm, cur_listing, cur_bytes); 
          buf += cur_bytes;
          tmap_refseq_destroy(refseq);
      } 
  }

  // pack the bwt 
  for(i=0;i<2;i++) { // forward/reverse
      cur_listing = (0 == i) ? TMAP_SHM_LISTING_BWT : TMAP_SHM_LISTING_REV_BWT;
      if(listing & cur_listing) {
          tmap_progress_print("packing %s bwt", (0 == i) ? "forward" : "reverse");
          bwt = tmap_bwt_read(fn_fasta, i);
          cur_bytes = tmap_bwt_shm_num_bytes(bwt);
          tmap_bwt_shm_pack(bwt, buf);
          tmap_shm_add_listing(shm, cur_listing, cur_bytes); 
          buf += cur_bytes;
          tmap_bwt_destroy(bwt);
      } 
  }

  // pack the SA
  for(i=0;i<2;i++) { // forward/reverse
      cur_listing = (0 == i) ? TMAP_SHM_LISTING_SA : TMAP_SHM_LISTING_REV_SA;
      if(listing & cur_listing) {
          tmap_progress_print("packing %s sa", (0 == i) ? "forward" : "reverse");
          sa = tmap_sa_read(fn_fasta, i);
          cur_bytes = tmap_sa_shm_num_bytes(sa);
          tmap_sa_shm_pack(sa, buf);
          tmap_shm_add_listing(shm, cur_listing, cur_bytes); 
          buf += cur_bytes;
          tmap_sa_destroy(sa);
      } 
  }

  tmap_progress_print2("shared memory packed");

  // set as ready
  tmap_shm_set_ready(shm);
  tmap_progress_print2("server started");

  // sleep while the stop signal has not been given
  while(TMAP_SHM_DEAD != tmap_shm_get_state(shm)) {
      sleep(TMAP_SERVER_SLEEP);
  }

  // destroy 
  tmap_progress_print("stopping server");
  tmap_shm_destroy(shm, 0);
  tmap_progress_print2("server stopped");
}

void
tmap_server_stop(key_t key)
{
  tmap_shm_t *shm = NULL;

  // get shared memory
  tmap_progress_print("retrieving shared memory");
  shm = tmap_shm_init(key, 0, 0);
  tmap_progress_print2("shared memory retrieved");

  // set as not ready
  tmap_progress_print("sending stop signal");
  tmap_shm_set_dead(shm);
  tmap_progress_print2("stop signal sent");

  // destroy 
  tmap_shm_destroy(shm, 0);
}

void
tmap_server_kill(key_t key)
{
  tmap_shm_t *shm = NULL;

  // get shared memory
  tmap_progress_print("retrieving shared memory");
  shm = tmap_shm_init(key, 0, 0);
  tmap_progress_print2("shared memory retrieved");

  // set as not ready
  tmap_progress_print("sending stop signal");
  tmap_shm_set_dead(shm);
  tmap_progress_print2("stop signal sent");

  // destroy 
  tmap_progress_print("stopping server");
  tmap_shm_destroy(shm, 1);
  tmap_progress_print2("server stopped");
}

static int
usage()
{
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Usage: %s server [options]", PACKAGE);
  tmap_file_fprintf(tmap_file_stderr, "\n");
  tmap_file_fprintf(tmap_file_stderr, "Options (optional):\n");
  tmap_file_fprintf(tmap_file_stderr, "         -f FILE     the FASTA reference file name\n");
  tmap_file_fprintf(tmap_file_stderr, "         -c STRING   server command [start|stop|kill]\n");
  tmap_file_fprintf(tmap_file_stderr, "         -k INT      the server key\n");
  tmap_file_fprintf(tmap_file_stderr, "         -a          load all\n");
  tmap_file_fprintf(tmap_file_stderr, "         -r          load the forward packed reference\n");
  tmap_file_fprintf(tmap_file_stderr, "         -R          load the reverse packed reference\n");
  tmap_file_fprintf(tmap_file_stderr, "         -b          load the forward bwt\n");
  tmap_file_fprintf(tmap_file_stderr, "         -B          load the reverse bwt\n");
  tmap_file_fprintf(tmap_file_stderr, "         -s          load the forward SA\n");
  tmap_file_fprintf(tmap_file_stderr, "         -S          load the reverse SA\n");
  tmap_file_fprintf(tmap_file_stderr, "         -v          print verbose progress information\n");
  tmap_file_fprintf(tmap_file_stderr, "         -h          print this message\n");
  tmap_file_fprintf(tmap_file_stderr, "\n");
  return 1;
}

int32_t
tmap_server_get_command_int(char *optarg)
{
  if(0 == strcmp("start", optarg)) return TMAP_SERVER_START;
  else if(0 == strcmp("stop", optarg)) return TMAP_SERVER_STOP;
  else if(0 == strcmp("kill", optarg)) return TMAP_SERVER_KILL;
  return TMAP_SERVER_UNKNOWN;
}

int
tmap_server_main(int argc, char *argv[])
{
  int c, help=0, cmd = -1;
  char *fn_fasta=NULL;
  key_t key=13;
  uint32_t listing = 0;

  while((c = getopt(argc, argv, "f:c:k:arRbBsSvh")) >= 0) {
      switch(c) {
        case 'a':
          listing |= TMAP_SHM_LISTING_REFSEQ;
          listing |= TMAP_SHM_LISTING_REV_REFSEQ;
          listing |= TMAP_SHM_LISTING_BWT;
          listing |= TMAP_SHM_LISTING_REV_BWT;
          listing |= TMAP_SHM_LISTING_SA;
          listing |= TMAP_SHM_LISTING_REV_SA;
          break;
        case 'f':
          fn_fasta = tmap_strdup(optarg); break;
        case 'c':
          cmd = tmap_server_get_command_int(optarg); break;
        case 'k': 
          key = atoi(optarg); break;
        case 'r':
          listing |= TMAP_SHM_LISTING_REFSEQ; break;
        case 'R':
          listing |= TMAP_SHM_LISTING_REV_REFSEQ; break;
        case 'b':
          listing |= TMAP_SHM_LISTING_BWT; break;
        case 'B':
          listing |= TMAP_SHM_LISTING_REV_BWT; break;
        case 's':
          listing |= TMAP_SHM_LISTING_SA; break;
        case 'S':
          listing |= TMAP_SHM_LISTING_REV_SA; break;
        case 'v': 
          tmap_progress_set_verbosity(1); break;
        case 'h': 
        default: 
          return usage();
      }
  }
  
  if(argc != optind || 1 == argc || 1 == help) {
      return usage();
  }

  if(TMAP_SERVER_START == cmd && 0 == listing) {
      tmap_error("no data structures to load", Exit, CommandLineArgument);
  }

  switch(cmd) {
    case TMAP_SERVER_START:
      tmap_server_start(fn_fasta, key, listing); break;
    case TMAP_SERVER_STOP:
      tmap_server_stop(key); break;
    case TMAP_SERVER_KILL:
      tmap_server_kill(key); break;
    case TMAP_SERVER_UNKNOWN:
    default:
      tmap_error("command not recognized", Exit, OutOfRange);
  }

  free(fn_fasta);

  return 0;
}
