#include <stdio.h>
#include <string.h>
#include <config.h>
#include <stdint.h>

#include "util/fmap_error.h"
#include "io/fmap_file.h"

extern int 
fmap_index(int argc, char *argv[]);
extern int 
fmap_refseq_fasta2pac_main(int argc, char *argv[]);
extern int 
fmap_bwt_pac2bwt_main(int argc, char *argv[]);
extern int
fmap_sa_bwt2sa_main(int argc, char *argv[]);
extern int 
fmap_map1(int argc, char *argv[]);
extern int 
fmap_debug_exact(int argc, char *argv[]);

fmap_file_t *fmap_file_stdout = NULL;
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
  fprintf(stderr, "Utilities:\n");
  fprintf(stderr, "         fasta2pac      creates the packed FASTA file\n");
  fprintf(stderr, "         pac2bwt        creates the BWT string file from the packed FASTA file\n");
  fprintf(stderr, "         bwt2sa         creates the SA file from the BWT string file\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Mapping:\n");
  fprintf(stderr, "         map1           mapping procedure #1\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Debugging:\n");
  fprintf(stderr, "         exact          perform simple exact matching\n");
  return 1;
}

int main(int argc, char *argv[])
{
  int ret = 0;
  
  fmap_file_stdout = fmap_file_fdopen(fileno(stdout), "w", FMAP_FILE_NO_COMPRESSION);
  fmap_file_stderr = fmap_file_fdopen(fileno(stderr), "w", FMAP_FILE_NO_COMPRESSION);

  if(argc < 2) ret = usage();
  else if (0 == strcmp("index", argv[1])) ret = fmap_index(argc-1, argv+1);
  else if (0 == strcmp("fasta2pac", argv[1])) ret = fmap_refseq_fasta2pac_main(argc-1, argv+1);
  else if (0 == strcmp("pac2bwt", argv[1])) ret = fmap_bwt_pac2bwt_main(argc-1, argv+1);
  else if (0 == strcmp("bwt2sa", argv[1])) ret = fmap_sa_bwt2sa_main(argc-1, argv+1);
  else if (0 == strcmp("map1", argv[1])) ret = fmap_map1(argc-1, argv+1);
  else if (0 == strcmp("exact", argv[1])) ret = fmap_debug_exact(argc-1, argv+1);
  else {
      fmap_error1(PACKAGE, "Unknown command", Exit, CommandLineArgument);
  }

  fmap_file_fclose(fmap_file_stdout);
  fmap_file_fclose(fmap_file_stderr);

  return ret;
}

