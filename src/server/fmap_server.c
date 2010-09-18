#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <unistd.h>

#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_progress.h"
#include "../index/fmap_refseq.h"
#include "fmap_shm.h"
#include "fmap_server.h"

static void
fmap_server_start(char *fn_fasta, key_t key)
{
  size_t n_bytes = 0;
  fmap_shm_t *shm = NULL;
  fmap_refseq_t *refseq = NULL;

  fmap_progress_print("starting server");

  // load
  fmap_progress_print("reading in data");
  refseq = fmap_refseq_read(fn_fasta, 0);
  fmap_progress_print2("data read in");

  // # of bytes
  n_bytes = fmap_refseq_shm_num_bytes(refseq);

  // get shared memory
  fmap_progress_print("retrieving shared memory");
  shm = fmap_shm_init(key, n_bytes, 1);
  fmap_progress_print2("shared memory retrieved");

  // store
  fmap_progress_print("packing shared memory");
  fmap_refseq_shm_pack(refseq, shm->buf);
  fmap_progress_print2("shared memory packed");

  // destroy
  fmap_refseq_destroy(refseq);

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

static void
fmap_server_test(key_t key)
{
  int32_t i;
  fmap_shm_t *shm = NULL;
  fmap_refseq_t *refseq = NULL;
  
  // get shared memory
  fmap_progress_print("retrieving shared memory");
  shm = fmap_shm_init(key, 0, 0);
  fmap_progress_print2("shared memory retrieved");

  fmap_progress_print("unpacking data");
  refseq = fmap_refseq_shm_unpack(shm->buf);
  fmap_progress_print2("data unpacked");

  fmap_file_fprintf(fmap_file_stderr, "num_annos=%d\n", refseq->num_annos);
  for(i=0;i<refseq->num_annos;i++) {
      fmap_file_fprintf(fmap_file_stderr, "[%d] [%s]\n", i, refseq->annos[i].name->s);
  }

  fmap_progress_print("destroying data");
  fmap_refseq_destroy(refseq);
  fmap_progress_print2("data destroyed");
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
  fmap_file_fprintf(fmap_file_stderr, "         -s          start the server\n");
  fmap_file_fprintf(fmap_file_stderr, "         -S          stop any server on this machine\n");
  fmap_file_fprintf(fmap_file_stderr, "         -z          remove the zombied shared memory segment\n");
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

  fmap_progress_set_start_time(clock());
  fmap_progress_set_command(argv[0]);

  while((c = getopt(argc, argv, "f:k:sSztvh")) >= 0) {
      switch(c) {
        case 'f':
          fn_fasta = fmap_strdup(optarg); break;
        case 'k': 
          key = atoi(optarg); break;
        case 's': 
          cmd = 0; break;
        case 'S': 
          cmd = 1; break;
        case 'z':
          cmd = 2; break;
        case 't':
          cmd = 3; break;
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

  switch(cmd) {
    case 0:
      fmap_server_start(fn_fasta, key); break;
    case 1:
      fmap_server_stop(key); break;
    case 2:
      fmap_server_kill(key); break;
    case 3:
      fmap_server_test(key); break;
    default:
      fmap_error("start not recognized", Exit, OutOfRange);
  }

  free(fn_fasta);

  return 0;
}
