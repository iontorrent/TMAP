#include <stdlib.h>
#include <config.h>
#ifdef HAVE_SAMTOOLS
#include <sam.h>
#endif
#include "../util/fmap_error.h"
#include "../util/fmap_alloc.h"
#include "../util/fmap_definitions.h"
#include "../index/fmap_refseq.h"
#include "../seq/fmap_sff.h"
#include "fmap_sfferr_aux.h"

#ifdef HAVE_SAMTOOLS

fmap_sfferr_aux_pr_flow_given_call_t *
fmap_sfferr_aux_pr_flow_given_call_init(fmap_sff_header_t *header)
{
  int32_t i, j;
  fmap_sfferr_aux_pr_flow_given_call_t *ptr = NULL;

  ptr = fmap_calloc(1, sizeof(fmap_sfferr_aux_pr_flow_given_call_t), "ptr");

  ptr->flow_length = header->flow_length;
  ptr->max_polymer = 12; // TODO: move this to a define block
  ptr->max_signal = 1 << 16;

  // allocate memory
  ptr->signal_counts = fmap_malloc(sizeof(uint16_t**)*ptr->flow_length, "ptr->signal_counts");
  for(i=0;i<ptr->flow_length;i++) {
      ptr->signal_counts[i] = fmap_malloc(sizeof(uint16_t*)*ptr->max_polymer, "ptr->signal_counts[i]");
      for(j=0;j<ptr->max_polymer;j++) {
          ptr->signal_counts[i][j] = fmap_calloc(ptr->max_signal, sizeof(uint16_t), "ptr->signal_counts[i][j]");
      }
  }

  return ptr;
} 

void
fmap_sfferr_aux_pr_flow_given_call_destroy(fmap_sfferr_aux_pr_flow_given_call_t *ptr)
{
  int32_t i, j;
  for(i=0;i<ptr->flow_length;i++) {
      for(j=0;j<ptr->max_polymer;j++) {
          free(ptr->signal_counts[i][j]);
      }
      free(ptr->signal_counts[i]);
  }
  free(ptr->signal_counts);
  free(ptr);
}

void
fmap_sfferr_aux_pr_flow_given_call_realloc(fmap_sfferr_aux_pr_flow_given_call_t *ptr, int32_t flow_length, int32_t max_polymer)
{
  int32_t i, j;

  if(ptr->flow_length < flow_length) {
      ptr->signal_counts = fmap_realloc(ptr->signal_counts, sizeof(uint16_t**)*ptr->flow_length, "ptr->signal_counts");
      for(i=ptr->flow_length;i<flow_length;i++) {
          ptr->signal_counts[i] = fmap_malloc(sizeof(uint16_t*)*ptr->max_polymer, "ptr->signal_counts[i]");
          for(j=0;j<ptr->max_polymer;j++) {
              ptr->signal_counts[i][j] = fmap_calloc(ptr->max_signal, sizeof(uint16_t), "ptr->signal_counts[i][j]");
          }
      }
      ptr->flow_length = flow_length;
  }
  if(ptr->max_polymer < max_polymer) {
      for(i=0;i<ptr->flow_length;i++) {
          ptr->signal_counts[i] = fmap_realloc(ptr->signal_counts[i], sizeof(uint16_t*)*max_polymer, "ptr->signal_counts[i]");
          for(j=ptr->max_polymer;j<max_polymer;j++) {
              ptr->signal_counts[i][j] = fmap_calloc(ptr->max_signal, sizeof(uint16_t), "ptr->signal_counts[i][j]");
          }
      }
      ptr->max_polymer = max_polymer;
  }
}

fmap_sfferr_aux_pr_flow_given_call_t *
fmap_sfferr_aux_pr_flow_given_call_read(fmap_file_t *fp)
{
  fmap_sfferr_aux_pr_flow_given_call_t *ptr;
  int32_t i, j;
  
  ptr = fmap_calloc(1, sizeof(fmap_sfferr_aux_pr_flow_given_call_t), "ptr");

  if(1 != fmap_file_fread(&ptr->flow_length, sizeof(int32_t), 1, fp)
     || 1 != fmap_file_fread(&ptr->max_polymer, sizeof(int32_t), 1, fp)
     || 1 != fmap_file_fread(&ptr->max_signal, sizeof(int32_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }

  ptr->signal_counts = fmap_malloc(sizeof(uint16_t**)*ptr->flow_length, "ptr->signal_counts");
  for(i=0;i<ptr->flow_length;i++) {
      ptr->signal_counts[i] = fmap_malloc(sizeof(uint16_t*)*ptr->max_polymer, "ptr->signal_counts[i]");
      for(j=0;j<ptr->max_polymer;j++) {
          ptr->signal_counts[i][j] = fmap_calloc(ptr->max_signal, sizeof(uint16_t), "ptr->signal_counts[i][j]");
          if(ptr->max_signal != fmap_file_fread(ptr->signal_counts[i][j], sizeof(uint16_t), ptr->max_signal, fp)) {
              fmap_error(NULL, Exit, WriteFileError);
          }
      }
  }

  return ptr;
} 

void
fmap_sfferr_aux_pr_flow_given_call_write(fmap_sfferr_aux_pr_flow_given_call_t *ptr, fmap_file_t *fp)
{
  int32_t i, j;

  if(1 != fmap_file_fwrite(&ptr->flow_length, sizeof(int32_t), 1, fp)
     || 1 != fmap_file_fwrite(&ptr->max_polymer, sizeof(int32_t), 1, fp)
     || 1 != fmap_file_fwrite(&ptr->max_signal, sizeof(int32_t), 1, fp)) {
      fmap_error(NULL, Exit, WriteFileError);
  }

  for(i=0;i<ptr->flow_length;i++) {
      for(j=0;j<ptr->max_polymer;j++) {
          if(ptr->max_signal != fmap_file_fwrite(ptr->signal_counts[i][j], sizeof(uint16_t), ptr->max_signal, fp)) {
              fmap_error(NULL, Exit, WriteFileError);
          }
      }
  }
} 

static void
fmap_sfferr_aux_pr_flow_given_call_add(fmap_sfferr_aux_pr_flow_given_call_t *ptr,
                                       fmap_sff_t *sff)
{
  int32_t i, j, is_int, n_bases = 0;

  // to int
  is_int = sff->is_int;
  fmap_sff_to_char(sff);

  for(i=j=0;i<sff->gheader->flow_length;i++) {
      n_bases = 0;
      while(j < sff->rheader->n_bases
            && sff->read->bases->s[j] == sff->gheader->flow->s[i]) {
          n_bases++;
          j++;
      }
      // i - flow index
      // j - base index
      // n_bases - the number bases called for this flow
      if(ptr->flow_length <= i || ptr->max_polymer <= n_bases) {
          // reallocate, if necessary
          fmap_sfferr_aux_pr_flow_given_call_realloc(ptr, i+1, n_bases+1);
      }
      ptr->signal_counts[i][n_bases][sff->read->flowgram[i]]++; // update the counts
  }

  // back to char if it was char
  if(1 == is_int) fmap_sff_to_int(sff);
}

void
fmap_sffer_aux_pr_call_given_polymer_add(fmap_refseq_t *refseq, fmap_sff_t *sff, bam1_t *bam)
{
  // TODO
}

void
fmap_sfferr_aux(fmap_refseq_t *refseq, fmap_sff_t *sff, bam1_t *bam,
                fmap_sfferr_aux_pr_flow_given_call_t *ptr1
                )
{
  // Pr(flow signal | polymer call, flow index)
  if(NULL != ptr1) fmap_sfferr_aux_pr_flow_given_call_add(ptr1, sff);

  // Pr(polymer call | reference polymer, flow index)
  // TODO
}

#endif /* HAVE_SAMTOOLS */
