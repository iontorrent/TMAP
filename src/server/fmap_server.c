#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <unistd.h>
#include <signal.h>

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
  fmap_error(NULL, Warn, SigInt); // warn
  if(NULL != fmap_server_shm_ptr) {
      // try to destroy the shared memory
      fmap_shm_set_dead(fmap_server_shm_ptr);
      fmap_shm_destroy(fmap_server_shm_ptr, 0);
      fmap_server_shm_ptr = NULL;
  }
  fmap_error(NULL, Exit, SigInt); //exit
}

static void
fmap_server_set_sigint(fmap_shm_t *shm)
{
  fmap_server_shm_ptr = shm;
  signal(SIGINT, fmap_server_sigint); 
}

// TODO: avoid loading twice
static void
fmap_server_start(char *fn_fasta, key_t key, uint32_t listing)
{
  int32_t i, j;
  fmap_shm_t *shm = NULL;
  fmap_refseq_t *refseq = NULL;
  fmap_bwt_t *bwt = NULL;
  fmap_sa_t *sa = NULL;
  uint8_t *buf = NULL;

  fmap_progress_print("starting server");

  // Two passes
  // - pass 1: get the number of bytes
  // - pass 2: store the memory
  for(i=0;i<2;i++) {
      size_t n_bytes = 0, cur_bytes = 0;
      if(0 == i) {
          fmap_progress_print("reading in data");
      }
      else {
          fmap_progress_print("packing shared memory");
          buf = (uint8_t*)shm->buf;
      }
      uint32_t cur_listing = 0;
      // reference sequence
      for(j=0;j<2;j++) { // forward/reverse
          cur_listing = (0 == j) ? FMAP_SERVER_LISTING_REFSEQ : FMAP_SERVER_LISTING_REV_REFSEQ;
          if(listing & cur_listing) {
              refseq = fmap_refseq_read(fn_fasta, j);
              cur_bytes = fmap_refseq_shm_num_bytes(refseq);
              if(1 == i) { // add to the shared memory
                  fmap_refseq_shm_pack(refseq, buf);
                  fmap_shm_add_listing(shm, cur_listing, cur_bytes); 
                  buf += cur_bytes;
              }
              fmap_refseq_destroy(refseq);
              n_bytes += cur_bytes;
          } 
      }
      // bwt 
      for(j=0;j<2;j++) { // forward/reverse
          cur_listing = (0 == j) ? FMAP_SERVER_LISTING_BWT : FMAP_SERVER_LISTING_REV_BWT;
          if(listing & cur_listing) {
              bwt = fmap_bwt_read(fn_fasta, j);
              cur_bytes = fmap_bwt_shm_num_bytes(bwt);
              if(1 == i) { // add to the shared memory
                  fmap_bwt_shm_pack(bwt, buf);
                  fmap_shm_add_listing(shm, cur_listing, cur_bytes); 
                  buf += cur_bytes;
              }
              fmap_bwt_destroy(bwt);
              n_bytes += cur_bytes;
          } 
      }
      // SA
      for(j=0;j<2;j++) { // forward/reverse
          cur_listing = (0 == j) ? FMAP_SERVER_LISTING_SA : FMAP_SERVER_LISTING_REV_SA;
          if(listing & cur_listing) {
              sa = fmap_sa_read(fn_fasta, j);
              cur_bytes = fmap_sa_shm_num_bytes(sa);
              if(1 == i) { // add to the shared memory
                  fmap_sa_shm_pack(sa, buf);
                  fmap_shm_add_listing(shm, cur_listing, cur_bytes); 
                  buf += cur_bytes;
              }
              n_bytes += cur_bytes;
              fmap_sa_destroy(sa);
          } 
      }
      if(0 == i) { // get the shared memory
          fmap_progress_print2("data read in");
          fmap_progress_print("retrieving shared memory [%llu bytes]", (long long unsigned int)n_bytes);
          shm = fmap_shm_init(key, n_bytes, 1);
          fmap_progress_print2("shared memory retrieved");
      
          // catch a ctrl-c signal
          fmap_server_set_sigint(shm);
      }
      else {
          fmap_progress_print2("shared memory packed");
      }
  } // i

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
  fmap_file_fprintf(fmap_file_stderr, "         -f FILE     the FASTA file name to index\n");
  fmap_file_fprintf(fmap_file_stderr, "         -k INT      the server key\n");
  fmap_file_fprintf(fmap_file_stderr, "         -0          stop the server\n");
  fmap_file_fprintf(fmap_file_stderr, "         -1          start the server\n");
  fmap_file_fprintf(fmap_file_stderr, "         -z          remove the zombied shared memory segment\n");
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

int
fmap_server_main(int argc, char *argv[])
{
  int c, help=0, cmd = 1;
  char *fn_fasta=NULL;
  key_t key=13;
  uint32_t listing = 0;

  fmap_progress_set_start_time(clock());
  fmap_progress_set_command(argv[0]);

  while((c = getopt(argc, argv, "f:k:01zrRbBsSvh")) >= 0) {
      switch(c) {
        case 'f':
          fn_fasta = fmap_strdup(optarg); break;
        case 'k': 
          key = atoi(optarg); break;
        case '0': 
          cmd = 0; break;
        case '1': 
          cmd = 1; break;
        case 'z':
          cmd = 2; break;
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

  if(1 == cmd && 0 == listing) {
      fmap_error("no data structures to load", Exit, CommandLineArgument);
  }

  switch(cmd) {
    case 0:
      fmap_server_stop(key); break;
    case 1:
      fmap_server_start(fn_fasta, key, listing); break;
    case 2:
      fmap_server_kill(key); break;
    default:
      fmap_error("start not recognized", Exit, OutOfRange);
  }

  free(fn_fasta);

  return 0;
}
