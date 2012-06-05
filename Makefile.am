bin_PROGRAMS = tmap

GLOBAL_SOURCES = \
				 src/util/tmap_alloc.h src/util/tmap_alloc.c \
				 src/util/tmap_definitions.h src/util/tmap_definitions.c \
				 src/util/tmap_error.h src/util/tmap_error.c \
				 src/util/tmap_rand.h src/util/tmap_rand.c \
				 src/util/tmap_time.h src/util/tmap_time.c \
				 src/util/tmap_fibheap.h src/util/tmap_fibheap.c \
				 src/util/tmap_hash.h \
				 src/util/tmap_progress.h src/util/tmap_progress.c \
				 src/util/tmap_string.h src/util/tmap_string.c \
				 src/util/tmap_levenshtein.h src/util/tmap_levenshtein.c \
				 src/util/tmap_sort.h \
				 src/util/tmap_vec.h \
				 src/util/tmap_sam_convert.h src/util/tmap_sam_convert.c \
				 src/seq/tmap_fq.h src/seq/tmap_fq.c \
				 src/seq/tmap_sff.h src/seq/tmap_sff.c \
				 src/seq/tmap_sam.h src/seq/tmap_sam.c \
				 src/seq/tmap_seq.h src/seq/tmap_seq.c \
				 src/seq/tmap_seqs.h src/seq/tmap_seqs.c \
				 src/io/tmap_file.h src/io/tmap_file.c \
				 src/io/tmap_fq_io.h src/io/tmap_fq_io.c \
				 src/io/tmap_sff_io.h src/io/tmap_sff_io.c \
				 src/io/tmap_sam_io.h src/io/tmap_sam_io.c \
				 src/io/tmap_seq_io.h src/io/tmap_seq_io.c \
				 src/io/tmap_seqs_io.h src/io/tmap_seqs_io.c \
				 src/index/tmap_bwt.h src/index/tmap_bwt.c \
				 src/index/tmap_bwt_aux.h src/index/tmap_bwt_aux.c \
				 src/index/tmap_bwtl.h src/index/tmap_bwtl.c \
				 src/index/tmap_bwt_gen.h src/index/tmap_bwt_gen.c \
				 src/index/tmap_bwt_match.h src/index/tmap_bwt_match.c \
				 src/index/tmap_bwt_match_hash.h src/index/tmap_bwt_match_hash.c \
				 src/index/tmap_bwt_smem.h src/index/tmap_bwt_smem.c \
				 src/index/tmap_index.h src/index/tmap_index.c \
				 src/index/tmap_refseq.h src/index/tmap_refseq.c \
				 src/index/tmap_sa.c src/index/tmap_sa.h \
				 src/index/tmap_sa_aux.c src/index/tmap_sa_aux.h \
				 src/index/tmap_index_size.h src/index/tmap_index_size.c \
				 src/sw/tmap_sw.h src/sw/tmap_sw.c \
				 src/sw/tmap_fsw.h src/sw/tmap_fsw.c \
				 src/sw/tmap_vsw_definitions.h src/sw/tmap_vsw_definitions.c \
				 src/sw/tmap_vsw.h src/sw/tmap_vsw.c \
				 src/sw/lib/vsw.cpp src/sw/lib/vsw.h \
				 src/sw/lib/vsw16.cpp src/sw/lib/vsw16.h \
				 src/sw/lib/sw-vector.cpp src/sw/lib/sw-vector.h \
				 src/sw/lib/Solution.cpp src/sw/lib/Solution.h \
				 src/sw/lib/Solution1.cpp src/sw/lib/Solution1.h \
				 src/sw/lib/Solution2.cpp src/sw/lib/Solution2.h \
				 src/sw/lib/Solution3.cpp src/sw/lib/Solution3.h \
				 src/sw/lib/Solution4.cpp src/sw/lib/Solution4.h \
				 src/sw/lib/Solution5.cpp src/sw/lib/Solution5.h \
				 src/sw/lib/Solution6.cpp src/sw/lib/Solution6.h \
				 src/sw/lib/Solution7.cpp src/sw/lib/Solution7.h \
				 src/sw/lib/Solution8.cpp src/sw/lib/Solution8.h \
				 src/sw/lib/Solution9.cpp src/sw/lib/Solution9.h \
				 src/sw/lib/Solution10.cpp src/sw/lib/Solution10.h \
				 src/sw/lib/AffineSWOptimization.cpp src/sw/lib/AffineSWOptimization.h \
				 src/sw/lib/AffineSWOptimizationHash.cpp src/sw/lib/AffineSWOptimizationHash.h \
				 src/sw/lib/AffineSWOptimizationWrapper.cpp src/sw/lib/AffineSWOptimizationWrapper.h \
				 src/sw/tmap_vsw_bm.c 

samtools_SOURCES = \
				   src/samtools/bam.c src/samtools/bam.h \
				   src/samtools/sam.c src/samtools/sam.h \
				   src/samtools/sam_header.c src/samtools/sam_header.h \
				   src/samtools/bam_import.c \
				   src/samtools/bam_aux.c \
				   src/samtools/bam_pileup.c \
				   src/samtools/razf.c src/samtools/razf.h \
				   src/samtools/bgzf.c src/samtools/bgzf.h \
				   src/samtools/faidx.c src/samtools/faidx.h \
				   src/samtools/kstring.c src/samtools/kstring.h

tmap_SOURCES = \
			   ${GLOBAL_SOURCES} \
			   ${samtools_SOURCES} \
			   src/server/tmap_server.h src/server/tmap_server.c \
			   src/server/tmap_shm.h src/server/tmap_shm.c \
			   src/sam2fs/tmap_sam2fs_aux.h src/sam2fs/tmap_sam2fs_aux.c \
			   src/sam2fs/tmap_sam2fs.h src/sam2fs/tmap_sam2fs.c \
			   src/index/tmap_index_speed.h src/index/tmap_index_speed.c \
			   src/index/tmap_bwt_check.c src/index/tmap_bwt_check.h \
			   src/index/tmap_bwt_compare.c src/index/tmap_bwt_compare.h \
			   src/map/util/tmap_map_opt.h src/map/util/tmap_map_opt.c \
			   src/map/util/tmap_map_stats.h src/map/util/tmap_map_stats.c \
			   src/map/util/tmap_map_util.h src/map/util/tmap_map_util.c \
			   src/map/pairing/tmap_map_pairing.h src/map/pairing/tmap_map_pairing.c \
			   src/map/map1/tmap_map1.h src/map/map1/tmap_map1.c \
			   src/map/map1/tmap_map1_aux.h src/map/map1/tmap_map1_aux.c \
			   src/map/map2/tmap_map2.h src/map/map2/tmap_map2.c \
			   src/map/map2/tmap_map2_aux.h src/map/map2/tmap_map2_aux.c \
			   src/map/map2/tmap_map2_chain.h src/map/map2/tmap_map2_chain.c \
			   src/map/map2/tmap_map2_core.h src/map/map2/tmap_map2_core.c \
			   src/map/map2/tmap_map2_mempool.h src/map/map2/tmap_map2_mempool.c \
			   src/map/map3/tmap_map3.h src/map/map3/tmap_map3.c \
			   src/map/map3/tmap_map3_aux.h src/map/map3/tmap_map3_aux.c \
			   src/map/map4/tmap_map4.h src/map/map4/tmap_map4.c \
			   src/map/map4/tmap_map4_aux.h src/map/map4/tmap_map4_aux.c \
			   src/map/mapvsw/tmap_map_vsw.h src/map/mapvsw/tmap_map_vsw.c \
			   src/map/mapall/tmap_map_all.h src/map/mapall/tmap_map_all.c \
			   src/map/tmap_map_driver.h src/map/tmap_map_driver.c \
			   src/tmap_main.h src/tmap_main.c

tmap_LDADD = 

#noinst_LIBRARIES = libtmap.a
#libtmap_a_SOURCES = ${GLOBAL_SOURCES} 

EXTRA_DIST = LICENSE \
			 autogen.sh \
			 tmap.git.rev \
			 doc/Makefile

SUBDIRS = . 

docdir = ${datadir}/doc/${PACKAGE}
dist_doc_DATA = LICENSE doc/tmap-book.pdf
