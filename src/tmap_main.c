/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <config.h>
#include <stdint.h>

#include "util/tmap_error.h"
#include "util/tmap_alloc.h"
#include "util/tmap_progress.h"
#include "util/tmap_levenshtein.h"
#include "util/tmap_definitions.h"
#include "io/tmap_file.h"
#include "samtools/bam.h"
#include "tmap_main.h"

tmap_file_t *tmap_file_stdout;
tmap_file_t *tmap_file_stderr;

static int 
tmap_usage(int argc, char *argv[]);

static tmap_command_t commands[] = {
      {tmap_index, "index", "creates the packed FASTA, BWT string, and SA files", TMAP_COMMAND_PREPROCESSING},
      {tmap_server_main, "server", "creates a mapping server", TMAP_COMMAND_SERVER},
      {tmap_map1_main, "map1", "mapping procedure #1 (bwa-short variant)", TMAP_COMMAND_MAPPING},
      {tmap_map2_main, "map2", "mapping procedure #2 (bwa-long/BWASW variant)", TMAP_COMMAND_MAPPING},
      {tmap_map3_main, "map3", "mapping procedure #3 (k-mer lookup)", TMAP_COMMAND_MAPPING},
      {tmap_map4_main, "map4", "mapping procedure #4 (bwa fastmap variant)", TMAP_COMMAND_MAPPING},
      {tmap_map_vsw_main, "mapvsw", "mapping procedure vectorized smith waterman", TMAP_COMMAND_MAPPING},
      {tmap_map_all_main, "mapall", "multi-mapping procedure", TMAP_COMMAND_MAPPING},
      {tmap_refseq_fasta2pac_main, "fasta2pac", "creates the packed FASTA file", TMAP_COMMAND_UTILITIES},
      {tmap_bwt_pac2bwt_main, "pac2bwt", "creates the BWT string file from the packed FASTA file", TMAP_COMMAND_UTILITIES},
      {tmap_sa_bwt2sa_main, "bwt2sa", "creates the SA file from the BWT string file", TMAP_COMMAND_UTILITIES},
      {tmap_seq_io_sff2fq_main, "sff2fq", "converts a SFF file to a FASTQ file", TMAP_COMMAND_UTILITIES},
      {tmap_seqs_io_sff2sam_main, "sff2sam", "converts a SFF file to a SAM file", TMAP_COMMAND_UTILITIES},
      {tmap_refseq_refinfo_main, "refinfo", "prints information about the reference", TMAP_COMMAND_UTILITIES},
      {tmap_refseq_pac2fasta_main, "pac2fasta", "converts a packed FASTA to a FASTA file", TMAP_COMMAND_UTILITIES},
      {tmap_bwt_bwtupdate_main, "bwtupdate", "updates the bwt hash width", TMAP_COMMAND_UTILITIES},
      {tmap_index_size, "indexsize", "gives the index size in bytes", TMAP_COMMAND_UTILITIES},
      {tmap_sam2fs_main, "sam2fs", "pretty print SAM records in flow space", TMAP_COMMAND_UTILITIES},
      {samtools_main, "samtools", "samtools (Tools for alignments in the SAM format)", TMAP_COMMAND_UTILITIES},
      {bcftools_main, "bcftools", "bcftools (Tools for data in the VCF/BCF formats)", TMAP_COMMAND_UTILITIES},
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
      {tmap_fsw_main, "fsw", "perform flow Smith-Waterman", TMAP_COMMAND_DEBUG},
      {tmap_index_speed, "indexspeed", "index performance benchmarks", TMAP_COMMAND_DEBUG},
      {tmap_bwt_check, "bwtcheck", "check the consistency of the BWT", TMAP_COMMAND_DEBUG},
      {tmap_bwt_compare, "bwtcompare", "compare two BWTs", TMAP_COMMAND_DEBUG},
      {tmap_vswbm_main, "vswbm", "VSW benchmarks", TMAP_COMMAND_DEBUG},
#endif 
      {tmap_version, "--version", "prints the TMAP version", TMAP_COMMAND_NONE},
      {tmap_version, "-v", "prints the TMAP version", TMAP_COMMAND_NONE},
      {tmap_usage, "--help", "prints this message", TMAP_COMMAND_NONE},
      {tmap_usage, "-h", "prints this message", TMAP_COMMAND_NONE},
      {NULL, NULL, NULL, -1}
};

static int
tmap_usage(int argc, char *argv[])
{
  tmap_command_t *c = NULL;
  int i, t;

  tmap_version(argc, argv);

  c = commands;
  t = -1;
  while(0 <= c->type) {
      if(c->type != t) {
          if(0 <= t && TMAP_COMMAND_NONE != c->type) fprintf(stderr, "\n");
          switch(c->type) {
            case TMAP_COMMAND_PREPROCESSING:
              fprintf(stderr, "%sPre-processing:%s\n", KRED, KNRM);
              break;
            case TMAP_COMMAND_SERVER:
              fprintf(stderr, "%sServer:%s\n", KRED, KNRM);
              break;
            case TMAP_COMMAND_MAPPING:
              fprintf(stderr, "%sMapping:%s\n", KRED, KNRM);
              break;
            case TMAP_COMMAND_UTILITIES:
              fprintf(stderr, "%sUtilities:%s\n", KRED, KNRM);
              break;
#ifdef ENABLE_TMAP_DEBUG_FUNCTIONS
            case TMAP_COMMAND_DEBUG:
              fprintf(stderr, "%sDebugging:%s\n", KRED, KNRM);
              break;
#endif
            case TMAP_COMMAND_NONE:
              break;
            default:
              fprintf(stderr, "c->type=%d\n", c->type);
              tmap_bug();
          }
          t = c->type;
      }
            
      if(c->type != TMAP_COMMAND_NONE) {
          fprintf(stderr, "         %s%s%s", KCYN, c->name, KNRM);
          for(i=strlen(c->name);i<16;i++) fputc(' ', stderr);
          fprintf(stderr, "%s%s%s\n", KWHT, c->help, KNRM);
      }
      c++;
  }
  return 1;
}

#define __get_distance(_a) ((int)(_a&0xffffffff))
#define __get_name_idx(_a) ((int)(_a>>32))

/* An empirically derived magic number */
#define TMAP_HELP_SIMILARITY_FLOOR 7
#define TMAP_HELP_SIMILAR_ENOUGH(x) ((x) < TMAP_HELP_SIMILARITY_FLOOR)

static int
tmap_prefixcmp(const char *str, const char *prefix)
{
  int32_t i;
  for(i=0;i<strlen(str) && i<strlen(prefix);i++) {
      if(str[i] != prefix[i])
        return (int)prefix[i] - (int)str[i];
  }
  //fprintf(stderr, "i=%d str=%s prefix=%s\n", i, str, prefix);
  return 0;
}

// NB: would require sorting to get the commands sorted
void
tmap_help_unknown_cmd(const char *cmd)
{
  int32_t i, n, best=INT32_MAX, best_n=0;
  uint64_t *distances = NULL;
  tmap_command_t *c = NULL;
  
  // get # of commands
  n = 0;
  c = commands;
  while(0 <= c->type) {
      n++;
      c++;
  }

  distances = tmap_malloc(sizeof(uint64_t)*n, "distances");

  for(i=0;0 <= commands[i].type;i++) {
      if(0 == strcmp(cmd, commands[i].name)) tmap_bug();
      if(!tmap_prefixcmp(commands[i].name, cmd)) {
          distances[i] = ((uint64_t)i << 32); // zero score
      }
      else {
          distances[i] = ((uint64_t)i << 32) | (tmap_levenshtein(cmd, commands[i].name, 0, 2, 1, 4) + 1); // pack
      }
      if(__get_distance(distances[i]) < best) {
          best = __get_distance(distances[i]);
          best_n = 0;
      }
      else if(__get_distance(distances[i]) == best) {
          best_n++;
      }
      //fprintf(stderr, "%s -> %u\n", commands[distances[i]>>32].name, (__get_distance(distances[i])));
  }

  if(0 == best && n == best_n) { // matches everything
      best = TMAP_HELP_SIMILARITY_FLOOR + 1; 
  }
  n = i;

  // output similar matches
  fprintf(stderr, "%s: '%s' is not a tmap command.  See 'tmap --help'.\n", PACKAGE, cmd);
  if(TMAP_HELP_SIMILAR_ENOUGH(best)) {
      fprintf(stderr, "\nDid you mean %s?\n",
              i < 2 ? "this": "one of these");
      for (i = 0; i < n; i++)
        if(best == __get_distance(distances[i])) {
            fprintf(stderr, "\t%s\n", commands[__get_name_idx(distances[i])].name);
        }
  }

  free(distances);
}

int 
main(int argc, char *argv[])
{
  int ret = 0;
  tmap_command_t *c = NULL;

  if(argc < 2) {
      return tmap_usage(argc, argv);
  }
  else {
      tmap_file_stderr = tmap_file_fdopen(fileno(stderr), "w", TMAP_FILE_NO_COMPRESSION); // set stderr
      tmap_progress_set_command(argv[1]); // set output progress
      tmap_progress_set_start_time(); // set start time
      bam_verbose = 0; // No verbosity in SAMtools

      c = commands;
      while(0 <= c->type) {
          if (0 == strcmp(c->name, argv[1])) {
              ret = c->func(argc-1, argv+1);
              break;
          }
          c++;
      }
      if(c->type < 0) {
          tmap_help_unknown_cmd(argv[1]);
          tmap_error1(PACKAGE, "Unknown command", Exit, CommandLineArgument);
      }

      // NB: do not close the underlying stderr stream!
      tmap_file_fclose1(tmap_file_stderr, 0);
  }

  return ret;
}
