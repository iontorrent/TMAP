#ifndef FMAP_MAIN_H_
#define FMAP_MAIN_H_

/*! 
  Main FMAP Function
  */

extern int 
fmap_index(int argc, char *argv[]);
extern int 
fmap_server_main(int argc, char *argv[]);
extern int 
fmap_map1_main(int argc, char *argv[]);
extern int
fmap_map2_main(int argc, char *argv[]);
extern int
fmap_map3_main(int argc, char *argv[]);
extern int 
fmap_refseq_fasta2pac_main(int argc, char *argv[]);
extern int 
fmap_bwt_pac2bwt_main(int argc, char *argv[]);
extern int
fmap_sa_bwt2sa_main(int argc, char *argv[]);
extern int
fmap_seq_io_sff2fq_main(int argc, char *argv[]);
#ifdef HAVE_SAMTOOLS
extern int
fmap_sam2fs_main(int argc, char *argv[]);
#endif
extern int 
fmap_debug_exact(int argc, char *argv[]);
extern int
fmap_fsw_main(int argc, char *argv[]);

#endif
