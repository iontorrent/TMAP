#ifndef FMAP_SEQ_H_
#define FMAP_SEQ_H_

#include "fmap_fq.h"
#include "fmap_sff.h"

/*! @enum  
  @constant  FMAP_SEQ_TYPE_SEQ  FASTA/FASTQ input/output
  @constant  FMAP_SEQ_TYPE_SFF  SFF input/output
 */
enum {
    FMAP_SEQ_TYPE_NOTYPE = -1,
    FMAP_SEQ_TYPE_FQ = 0,
    FMAP_SEQ_TYPE_SFF = 1
};

/*! @typedef 
  @field  type  the type associated with this structure
  @field  data  pointer to the particular read data structure
*/
typedef struct {
    int8_t type;
    union {
        fmap_fq_t *fq;
        fmap_sff_t *sff;
    } data;
} fmap_seq_t;

/*! @function
  @abstract
  @param  type  the type associated with this structure
  @return       pointer to the initialized memory 
*/
fmap_seq_t *
fmap_seq_init(int8_t type);

/*! @function
  @abstract
  @param  seq  pointer to the structure
  */
void
fmap_seq_destroy(fmap_seq_t *seq);

#endif
