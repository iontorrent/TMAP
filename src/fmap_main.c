#include <stdio.h>
#include <string.h>
#include <config.h>
#include <stdint.h>

#include "util/fmap_error.h"
#include "util/fmap_progress.h"
#include "io/fmap_file.h"
#include "fmap_main.h"

fmap_file_t *fmap_file_stdout = NULL; // do not initialize as this may used for output
fmap_file_t *fmap_file_stderr = NULL;

static int usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "%s:   flow mapper\n", PACKAGE);
#ifdef GIT_REV
  fprintf(stderr, "Version: %s git:%s\n", PACKAGE_VERSION, GIT_REV);
#else
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
#endif
  fprintf(stderr, "Contact: %s\n\n", PACKAGE_BUGREPORT);
  fprintf(stderr, "Usage:   %s <command> [options]\n\n", PACKAGE); 
  fprintf(stderr, "Pre-processing:\n");
  fprintf(stderr, "         index          creates the packed FASTA, BWT string, and SA files\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Server:\n");
  fprintf(stderr, "         server         creates a mapping server\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Mapping:\n");
  fprintf(stderr, "         map1           mapping procedure #1\n");
  fprintf(stderr, "         map2           mapping procedure #2\n");
  fprintf(stderr, "         map3           mapping procedure #3\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Utilities:\n");
  fprintf(stderr, "         fasta2pac      creates the packed FASTA file\n");
  fprintf(stderr, "         pac2bwt        creates the BWT string file from the packed FASTA file\n");
  fprintf(stderr, "         bwt2sa         creates the SA file from the BWT string file\n");
  fprintf(stderr, "         sff2fq         converts a SFF file to a FASTQ file\n");
#ifdef HAVE_SAMTOOLS
  fprintf(stderr, "         sfferr         creates an error profile frm an sff file and sam file\n");
#endif
  fprintf(stderr, "\n");
  fprintf(stderr, "Debugging:\n");
  fprintf(stderr, "         exact          perform simple exact matching\n");
  fprintf(stderr, "         fsw            perform flow Smith-Waterman\n");
  return 1;
}

int main(int argc, char *argv[])
{
  int ret = 0;

  if(argc < 2) {
      return usage();
  }
  else {

      fmap_file_stderr = fmap_file_fdopen(fileno(stderr), "w", FMAP_FILE_NO_COMPRESSION); // set stderr
      fmap_progress_set_command(argv[1]); // set output progress
      fmap_progress_set_start_time(clock()); // set start time

      if (0 == strcmp("index", argv[1])) ret = fmap_index(argc-1, argv+1);
      else if (0 == strcmp("server", argv[1])) ret = fmap_server_main(argc-1, argv+1);
      else if (0 == strcmp("map1", argv[1])) ret = fmap_map1_main(argc-1, argv+1);
      else if (0 == strcmp("map2", argv[1])) ret = fmap_map2_main(argc-1, argv+1);
      else if (0 == strcmp("map3", argv[1])) ret = fmap_map3_main(argc-1, argv+1);
      else if (0 == strcmp("fasta2pac", argv[1])) ret = fmap_refseq_fasta2pac_main(argc-1, argv+1);
      else if (0 == strcmp("pac2bwt", argv[1])) ret = fmap_bwt_pac2bwt_main(argc-1, argv+1);
      else if (0 == strcmp("bwt2sa", argv[1])) ret = fmap_sa_bwt2sa_main(argc-1, argv+1);
      else if (0 == strcmp("sff2fq", argv[1])) ret = fmap_seq_io_sff2fq_main(argc-1, argv+1);
#ifdef HAVE_SAMTOOLS
      else if (0 == strcmp("sfferr", argv[1])) ret = fmap_sfferr_main(argc-1, argv+1);
#endif
      else if (0 == strcmp("exact", argv[1])) ret = fmap_debug_exact(argc-1, argv+1);
      else if (0 == strcmp("fsw", argv[1])) ret = fmap_fsw_main(argc-1, argv+1);
      else {
          fmap_error1(PACKAGE, "Unknown command", Exit, CommandLineArgument);
      }

      fmap_file_fclose(fmap_file_stderr);
  }

  return ret;
}

