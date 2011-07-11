bin_PROGRAMS = tmap

GLOBAL_SOURCES = \
				 src/util/tmap_alloc.h src/util/tmap_alloc.c \
				 src/util/tmap_definitions.h src/util/tmap_definitions.c \
				 src/util/tmap_error.h src/util/tmap_error.c \
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
				 src/index/tmap_bwtl.h src/index/tmap_bwtl.c \
				 src/index/tmap_bwt_gen.h src/index/tmap_bwt_gen.c \
				 src/index/tmap_bwt_match.h src/index/tmap_bwt_match.c \
				 src/index/tmap_index.h src/index/tmap_index.c \
				 src/index/tmap_refseq.h src/index/tmap_refseq.c \
				 src/index/tmap_sa.c src/index/tmap_sa.h \
				 src/sw/tmap_sw.h src/sw/tmap_sw.c \
				 src/sw/tmap_fsw.h src/sw/tmap_fsw.c \
				 src/sw/tmap_vsw.h src/sw/tmap_vsw.c \
				 src/sw/tmap_vsw16.h src/sw/tmap_vsw16.c 

tmap_SOURCES = \
			   ${GLOBAL_SOURCES} \
			   src/server/tmap_server.h src/server/tmap_server.c \
			   src/server/tmap_shm.h src/server/tmap_shm.c \
			   src/sam2fs/tmap_sam2fs_aux.h src/sam2fs/tmap_sam2fs_aux.c \
			   src/sam2fs/tmap_sam2fs.h src/sam2fs/tmap_sam2fs.c \
			   src/debug/tmap_debug_exact.h src/debug/tmap_debug_exact.c \
			   src/map/tmap_map_vsw.h src/map/tmap_map_vsw.c \
			   src/map/tmap_map_util.h src/map/tmap_map_util.c \
			   src/map/tmap_map1.h src/map/tmap_map1.c \
			   src/map/tmap_map1_aux.h src/map/tmap_map1_aux.c \
			   src/map/tmap_map2.h src/map/tmap_map2.c \
			   src/map/tmap_map2_aux.h src/map/tmap_map2_aux.c \
			   src/map/tmap_map2_chain.h src/map/tmap_map2_chain.c \
			   src/map/tmap_map2_core.h src/map/tmap_map2_core.c \
			   src/map/tmap_map2_mempool.h src/map/tmap_map2_mempool.c \
			   src/map/tmap_map3.h src/map/tmap_map3.c \
			   src/map/tmap_map3_aux.h src/map/tmap_map3_aux.c \
			   src/map/tmap_map_all.h src/map/tmap_map_all.c \
			   src/map/tmap_map_driver.h src/map/tmap_map_driver.c \
			   src/tmap_main.h src/tmap_main.c

tmap_LDADD = 

#noinst_LIBRARIES = libtmap.a
#libtmap_a_SOURCES = ${GLOBAL_SOURCES} 

EXTRA_DIST = LICENSE \
			 autogen.sh 

SUBDIRS = .

docdir = ${datadir}/doc/${PACKAGE}
dist_doc_DATA = LICENSE doc/tmap-book.pdf
