bin_PROGRAMS = fmap

fmap_SOURCES = \
			   src/util/fmap_error.h src/util/fmap_error.c \
			   src/util/fmap_alloc.h src/util/fmap_alloc.c \
			   src/util/fmap_fibheap.h src/util/fmap_fibheap.c \
			   src/io/fmap_io.h src/io/fmap_io.c \
			   src/util/fmap_seq.h src/util/fmap_seq.c \
			   src/io/fmap_seq_io.h src/io/fmap_seq_io.c \
			   src/index/fmap_refseq.h src/index/fmap_refseq.c \
			   src/index/fmap_bwt.h src/index/fmap_bwt.c \
			   src/index/fmap_bwt_match.h src/index/fmap_bwt_match.c \
			   src/index/fmap_bwt_gen.h src/index/fmap_bwt_gen.c \
			   src/index/fmap_sa.c src/index/fmap_sa.h \
			   src/index/fmap_index.h src/index/fmap_index.c \
			   src/debug/fmap_debug_exact.h src/debug/fmap_debug_exact.c \
			   src/util/fmap_definitions.h src/util/fmap_definitions.c \
			   src/util/fmap_progress.h src/util/fmap_progress.c \
			   src/search/fmap_map1_aux.h src/search/fmap_map1_aux.c \
			   src/search/fmap_map1.h src/search/fmap_map1.c \
			   src/main.c

fmap_LDADD =

EXTRA_DIST = LICENSE \
			 autogen.sh 

SUBDIRS = .

docdir = ${datadir}/doc/${PACKAGE}
dist_doc_DATA = LICENSE 
