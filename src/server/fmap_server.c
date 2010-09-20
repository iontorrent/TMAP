#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <unistd.h>
#include <signal.h>
#include <string.h>

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_progress.h"
#include "../index/fmap_refseq.h"
#include "../index/fmap_bwt.h"
#include "../index/fmap_sa.h"
#include "fmap_shm.h"
#include "fmap_server.h"

static fmap_shm_t *fmap_server_shm_ptr = NULL;

static void
fmap_server_sigint(int signal)
{
  if(NULL != fmap_server_shm_ptr) {
      fmap_error("will try to destroy the shared memory", Warn, SigInt); // warn
      // try to destroy the shared memory
      fmap_shm_set_dead(fmap_server_shm_ptr);
      fmap_shm_destroy(fmap_server_shm_ptr, 0);
      fmap_server_shm_ptr = NULL;
      fmap_error("shared memory destroyed", Exit, SigInt); //exit
  }
  else {
      fmap_error(NULL, Exit, SigInt); //exit
  }
}

static void
fmap_server_set_sigint(fmap_shm_t *shm)
{
  fmap_server_shm_ptr = shm;
  signal(SIGINT, fmap_server_sigint); 
}

static void
fmap_server_start(char *fn_fasta, key_t key, uint32_t listing)
{
  int32_t i;
  fmap_shm_t *shm = NULL;
  fmap_refseq_t *refseq = NULL;
  fmap_bwt_t *bwt = NULL;
  fmap_sa_t *sa = NULL;
  uint8_t *buf = NULL;
  size_t n_bytes = 0, cur_bytes = 0;
  uint32_t cur_listing = 0;

  fmap_progress_print("starting server");

  // get data size
  n_bytes = 0;
  if(listing & FMAP_SERVER_LISTING_REFSEQ) {
      n_bytes += fmap_refseq_shm_read_num_bytes(fn_fasta, 0);
  }
  if(listing & FMAP_SERVER_LISTING_REV_REFSEQ) {
      n_bytes += fmap_refseq_shm_read_num_bytes(fn_fasta, 1);
  }
  if(listing & FMAP_SERVER_LISTING_BWT) {
      n_bytes += fmap_bwt_shm_read_num_bytes(fn_fasta, 0);
  }
  if(listing & FMAP_SERVER_LISTING_REV_BWT) {
      n_bytes += fmap_bwt_shm_read_num_bytes(fn_fasta, 1);
  }
  if(listing & FMAP_SERVER_LISTING_SA) {
      n_bytes += fmap_sa_shm_read_num_bytes(fn_fasta, 0);
  }
  if(listing & FMAP_SERVER_LISTING_REV_SA) {
      n_bytes += fmap_sa_shm_read_num_bytes(fn_fasta, 1);
  }

  // get shared memory
  fmap_progress_print("retrieving shared memory [%llu bytes]", (long long unsigned int)n_bytes);
  shm = fmap_shm_init(key, n_bytes, 1);
  fmap_progress_print2("shared memory retrieved");

  // catch a ctrl-c signal from now on
  fmap_server_set_sigint(shm);

  // pack the shared memory
  fmap_progress_print("packing shared memory");
  buf = (uint8_t*)shm->buf;

  // pack the reference sequence
  for(i=0;i<2;i++) { // forward/reverse
      cur_listing = (0 == i) ? FMAP_SERVER_LISTING_REFSEQ : FMAP_SERVER_LISTING_REV_REFSEQ;
      if(listing & cur_listing) {
          fmap_progress_print("packing %s reference", (0 == i) ? "forward" : "reverse");
          refseq = fmap_refseq_read(fn_fasta, i);
          cur_bytes = fmap_refseq_shm_num_bytes(refseq);
          fmap_refseq_shm_pack(refseq, buf);
          fmap_shm_add_listing(shm, cur_listing, cur_bytes); 
          buf += cur_bytes;
          fmap_refseq_destroy(refseq);
      } 
  }

  // pack the bwt 
  for(i=0;i<2;i++) { // forward/reverse
      cur_listing = (0 == i) ? FMAP_SERVER_LISTING_BWT : FMAP_SERVER_LISTING_REV_BWT;
      if(listing & cur_listing) {
          fmap_progress_print("packing %s bwt", (0 == i) ? "forward" : "reverse");
          bwt = fmap_bwt_read(fn_fasta, i);
          cur_bytes = fmap_bwt_shm_num_bytes(bwt);
          fmap_bwt_shm_pack(bwt, buf);
          fmap_shm_add_listing(shm, cur_listing, cur_bytes); 
          buf += cur_bytes;
          fmap_bwt_destroy(bwt);
      } 
  }

  // pack the SA
  for(i=0;i<2;i++) { // forward/reverse
      cur_listing = (0 == i) ? FMAP_SERVER_LISTING_SA : FMAP_SERVER_LISTING_REV_SA;
      if(listing & cur_listing) {
          fmap_progress_print("packing %s sa", (0 == i) ? "forward" : "reverse");
          sa = fmap_sa_read(fn_fasta, i);
          cur_bytes = fmap_sa_shm_num_bytes(sa);
          fmap_sa_shm_pack(sa, buf);
          fmap_shm_add_listing(shm, cur_listing, cur_bytes); 
          buf += cur_bytes;
          fmap_sa_destroy(sa);
      } 
  }

  fmap_progress_print2("shared memory packed");

  // set as ready
  fmap_shm_set_ready(shm);
  fmap_progress_print2("server started");

  // sleep while the stop signal has not been given
  while(FMAP_SHM_DEAD != fmap_shm_get_state(shm)) {
      sleep(FMAP_SERVER_SLEEP);
  }

  // destroy 
  fmap_progress_print("stopping server");
  fmap_shm_destroy(shm, 0);
  fmap_progress_print2("server stopped");
}

static void
fmap_server_stop(key_t key)
{
  fmap_shm_t *shm = NULL;

  // get shared memory
  fmap_progress_print("retrieving shared memory");
  shm = fmap_shm_init(key, 0, 0);
  fmap_progress_print2("shared memory retrieved");

  // set as not ready
  fmap_progress_print("sending stop signal");
  fmap_shm_set_dead(shm);
  fmap_progress_print2("stop signal sent");

  // destroy 
  fmap_shm_destroy(shm, 0);
}

static void
fmap_server_kill(key_t key)
{
  fmap_shm_t *shm = NULL;

  // get shared memory
  fmap_progress_print("retrieving shared memory");
  shm = fmap_shm_init(key, 0, 0);
  fmap_progress_print2("shared memory retrieved");

  // set as not ready
  fmap_progress_print("sending stop signal");
  fmap_shm_set_dead(shm);
  fmap_progress_print2("stop signal sent");

  // destroy 
  fmap_progress_print("stopping server");
  fmap_shm_destroy(shm, 1);
  fmap_progress_print2("server stopped");
}

static int
usage()
{
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Usage: %s server [options]", PACKAGE);
  fmap_file_fprintf(fmap_file_stderr, "\n");
  fmap_file_fprintf(fmap_file_stderr, "Options (optional):\n");
  fmap_file_fprintf(fmap_file_stderr, "         -f FILE     the FASTA reference file name\n");
  fmap_file_fprintf(fmap_file_stderr, "         -c STRING   server command [start|stop|kill]\n");
  fmap_file_fprintf(fmap_file_stderr, "         -k INT      the server key\n");
  fmap_file_fprintf(fmap_file_stderr, "         -r          load the forward packed reference\n");
  fmap_file_fprintf(fmap_file_stderr, "         -R          load the reverse packed reference\n");
  fmap_file_fprintf(fmap_file_stderr, "         -b          load the forward bwt\n");
  fmap_file_fprintf(fmap_file_stderr, "         -B          load the reverse bwt\n");
  fmap_file_fprintf(fmap_file_stderr, "         -s          load the forward SA\n");
  fmap_file_fprintf(fmap_file_stderr, "         -S          load the reverse SA\n");
  fmap_file_fprintf(fmap_file_stderr, "         -v          print verbose progress information\n");
  fmap_file_fprintf(fmap_file_stderr, "         -h          print this message\n");
  fmap_file_fprintf(fmap_file_stderr, "\n");
  return 1;
}

static int32_t
fmap_server_get_command_int(char *optarg)
{
  if(0 == strcmp("start", optarg)) return FMAP_SERVER_START;
  else if(0 == strcmp("stop", optarg)) return FMAP_SERVER_STOP;
  else if(0 == strcmp("kill", optarg)) return FMAP_SERVER_KILL;
  return FMAP_SERVER_UNKNOWN;
}

int
fmap_server_main(int argc, char *argv[])
{
  int c, help=0, cmd = -1;
  char *fn_fasta=NULL;
  key_t key=13;
  uint32_t listing = 0;

  fmap_progress_set_start_time(clock());
  fmap_progress_set_command(argv[0]);

  while((c = getopt(argc, argv, "f:a:k:rRbBsSvh")) >= 0) {
      switch(c) {
        case 'f':
          fn_fasta = fmap_strdup(optarg); break;
        case 'a':
          cmd = fmap_server_get_command_int(optarg); break;
        case 'k': 
          key = atoi(optarg); break;
        case 'r':
          listing |= FMAP_SERVER_LISTING_REFSEQ; break;
        case 'R':
          listing |= FMAP_SERVER_LISTING_REV_REFSEQ; break;
        case 'b':
          listing |= FMAP_SERVER_LISTING_BWT; break;
        case 'B':
          listing |= FMAP_SERVER_LISTING_REV_BWT; break;
        case 's':
          listing |= FMAP_SERVER_LISTING_SA; break;
        case 'S':
          listing |= FMAP_SERVER_LISTING_REV_SA; break;
        case 'v': 
          fmap_progress_set_verbosity(1); break;
        case 'h': 
        default: 
          return usage();
      }
  }
  if(argc != optind || 1 == argc || 1 == help) {
      return usage();
  }

  if(FMAP_SERVER_START == cmd && 0 == listing) {
      fmap_error("no data structures to load", Exit, CommandLineArgument);
  }

  switch(cmd) {
    case FMAP_SERVER_START:
      fmap_server_start(fn_fasta, key, listing); break;
    case FMAP_SERVER_STOP:
      fmap_server_stop(key); break;
    case FMAP_SERVER_KILL:
      fmap_server_kill(key); break;
    case FMAP_SERVER_UNKNOWN:
    default:
      fmap_error("command not recognized", Exit, OutOfRange);
  }

  free(fn_fasta);

  return 0;
}
