#include <stdio.h>
#include <string.h>
#include <config.h>
#include <stdint.h>

#include "util/tmap_error.h"
#include "util/tmap_progress.h"
#include "io/tmap_file.h"
#include "tmap_main.h"

tmap_file_t *tmap_file_stdout = NULL; // do not initialize as this may used for output
tmap_file_t *tmap_file_stderr = NULL;

static int version()
{
  fprintf(stdout, "\n");
  fprintf(stdout, "%s:   torrent mapper\n", PACKAGE);
#ifdef GIT_REV
  fprintf(stdout, "Version: %s git:%s\n", PACKAGE_VERSION, GIT_REV);
#else
  fprintf(stdout, "Version: %s\n", PACKAGE_VERSION);
#endif
  fprintf(stdout, "Contact: %s\n\n", PACKAGE_BUGREPORT);
  return 0;
}

static int usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "%s:   torrent mapper\n", PACKAGE);
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
  fprintf(stderr, "         mapall         multi-mapping procedure\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Utilities:\n");
  fprintf(stderr, "         fasta2pac      creates the packed FASTA file\n");
  fprintf(stderr, "         pac2bwt        creates the BWT string file from the packed FASTA file\n");
  fprintf(stderr, "         bwt2sa         creates the SA file from the BWT string file\n");
  fprintf(stderr, "         sff2fq         converts a SFF file to a FASTQ file\n");
#ifdef HAVE_SAMTOOLS
  fprintf(stderr, "         sam2fs         pretty print SAM records in flow space\n");
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

      tmap_file_stderr = tmap_file_fdopen(fileno(stderr), "w", TMAP_FILE_NO_COMPRESSION); // set stderr
      tmap_progress_set_command(argv[1]); // set output progress
      tmap_progress_set_start_time(clock()); // set start time

      if (0 == strcmp("index", argv[1])) ret = tmap_index(argc-1, argv+1);
      else if (0 == strcmp("server", argv[1])) ret = tmap_server_main(argc-1, argv+1);
      else if (0 == strcmp("map1", argv[1])) ret = tmap_map1_main(argc-1, argv+1);
      else if (0 == strcmp("map2", argv[1])) ret = tmap_map2_main(argc-1, argv+1);
      else if (0 == strcmp("map3", argv[1])) ret = tmap_map3_main(argc-1, argv+1);
      else if (0 == strcmp("mapall", argv[1])) ret = tmap_map_all_main(argc-1, argv+1);
      else if (0 == strcmp("fasta2pac", argv[1])) ret = tmap_refseq_fasta2pac_main(argc-1, argv+1);
      else if (0 == strcmp("pac2bwt", argv[1])) ret = tmap_bwt_pac2bwt_main(argc-1, argv+1);
      else if (0 == strcmp("bwt2sa", argv[1])) ret = tmap_sa_bwt2sa_main(argc-1, argv+1);
      else if (0 == strcmp("sff2fq", argv[1])) ret = tmap_seq_io_sff2fq_main(argc-1, argv+1);
#ifdef HAVE_SAMTOOLS
      else if (0 == strcmp("sam2fs", argv[1])) ret = tmap_sam2fs_main(argc-1, argv+1);
#endif
      else if (0 == strcmp("exact", argv[1])) ret = tmap_debug_exact(argc-1, argv+1);
      else if (0 == strcmp("fsw", argv[1])) ret = tmap_fsw_main(argc-1, argv+1);
      else if (0 == strcmp("--version", argv[1]) || 0 == strcmp("-v", argv[1])) ret = version();
      else if (0 == strcmp("--help", argv[1]) || 0 == strcmp("-h", argv[1])) ret = usage();
      else {
          tmap_error1(PACKAGE, "Unknown command", Exit, CommandLineArgument);
      }

      tmap_file_fclose(tmap_file_stderr);
  }

  return ret;
}

