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
				 src/util/tmap_sort.h \
				 src/util/tmap_vec.h \
				 src/util/tmap_sam_print.h src/util/tmap_sam_print.c \
				 src/seq/tmap_fq.h src/seq/tmap_fq.c \
				 src/seq/tmap_sff.h src/seq/tmap_sff.c \
				 src/seq/tmap_sam.h src/seq/tmap_sam.c \
				 src/seq/tmap_seq.h src/seq/tmap_seq.c \
				 src/io/tmap_file.h src/io/tmap_file.c \
				 src/io/tmap_fq_io.h src/io/tmap_fq_io.c \
				 src/io/tmap_sff_io.h src/io/tmap_sff_io.c \
				 src/io/tmap_sam_io.h src/io/tmap_sam_io.c \
				 src/io/tmap_seq_io.h src/io/tmap_seq_io.c \
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
				 src/index/tmap_index_size.h src/index/tmap_index_size.c \
				 src/sw/tmap_sw.h src/sw/tmap_sw.c \
				 src/sw/tmap_fsw.h src/sw/tmap_fsw.c \
				 src/sw/tmap_vsw.h src/sw/tmap_vsw.c \
				 src/sw/tmap_vsw16.h src/sw/tmap_vsw16.c \
				 src/sw/tmap_vsw_bm.c 

tmap_SOURCES = \
			   ${GLOBAL_SOURCES} \
			   src/server/tmap_server.h src/server/tmap_server.c \
			   src/server/tmap_shm.h src/server/tmap_shm.c \
			   src/sam2fs/tmap_sam2fs_aux.h src/sam2fs/tmap_sam2fs_aux.c \
			   src/sam2fs/tmap_sam2fs.h src/sam2fs/tmap_sam2fs.c \
			   src/debug/tmap_debug_exact.h src/debug/tmap_debug_exact.c \
			   src/index/tmap_index_speed.h src/index/tmap_index_speed.c \
			   src/index/tmap_bwt_check.c src/index/tmap_bwt_check.h \
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
			 tmap.git.rev

SUBDIRS = .

docdir = ${datadir}/doc/${PACKAGE}
dist_doc_DATA = LICENSE doc/tmap-book.pdf
