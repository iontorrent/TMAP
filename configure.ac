#  See LICENSE for details
#
#  The version number in AC_INIT is for the tmap package.
#  Developers - the version number should *NEVER* be changed except by, 
#  or with explicit permission from, Nils Homer.

AC_PREREQ(2.59);
AC_INIT([tmap],[0.0.12],[nils.homer@lifetech.com],[tmap])
AC_COPYRIGHT([See LICENSE for copyright information.])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([dist-bzip2 subdir-objects])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

# Sets build_cpu, build_vendor and build_os variables
AC_CANONICAL_BUILD

# Check for C and C++ compilers
AC_PROG_CC
AM_PROG_CC_C_O

AC_PROG_INSTALL
AC_GNU_SOURCE
#AC_PROG_RANLIB

# set CFLAGS and CXXFLAGS
default_CFLAGS="-Wall -g -O3 -pthread";# "-Wall -g -O2 -pthread";
extended_CFLAGS="";# "-m64 -D_FILE_OFFSET_BITS=64";

AC_MSG_CHECKING(whether this is a git repository)
gitrev="`git log 2>/dev/null | grep -m 1 commit | awk '{print $2}'`";
if test -n "${gitrev}"; then
    GITREV="${gitrev}";
    AC_REVISION(["${GITREV}"])
    AC_DEFINE_UNQUOTED(GIT_REV, "${GITREV}", [Define to be the git revision number])
    AC_MSG_RESULT(yes)
else
    AC_MSG_RESULT(no)
fi

AC_ARG_ENABLE(intel64, [  --enable-intel64        optimize for Intel64 CPU such as Xeon and Core2], [extended_CFLAGS="${extended_CFLAGS} -mtune=nocona"], [])
AC_ARG_ENABLE(bz2, [  --disable-bz2         use this option to disable bz2 support], [AC_DEFINE(DISABLE_BZ2,1,[Define 1 if we want to disable bz2 support])], [AC_CHECK_LIB([bz2], [BZ2_bzRead])])

CFLAGS="${default_CFLAGS} ${extended_CFLAGS}";

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
# Checks for library functions.
AC_FUNC_MALLOC 
AC_FUNC_REALLOC
AC_CHECK_LIB([z], [gzread])
AC_CHECK_LIB([bz2], [BZ2_bzRead])
AC_CHECK_LIB([m], [pow])
AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_FUNCS([pow strdup memset strchr strdup strstr memmove])

# Check types
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_SIZE_T

# Set proprocessor SIZEOF_short_int, SIZEOF_int, SIZEOF_long_int
AC_CHECK_SIZEOF([short int])
AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([long int])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h stdio.h stdint.h string.h sys/time.h unistd.h float.h zlib.h bzlib.h pthread.h limits.h netinet/in.h])

AC_C_INLINE

# test for SAMtools
AC_MSG_CHECKING(for samtools)
if test -d "samtools" ; then
	AC_MSG_RESULT(yes)
	LDFLAGS="${LDFLAGS} -L./samtools"
	LIBS="-lbam ${LIBS}";
	CFLAGS="${CFLAGS} -I./samtools";
	AC_DEFINE(HAVE_SAMTOOLS, 1, [Define 1 if the samtools library is present]) 
else 
	AC_MSG_RESULT(no)
	AC_MSG_WARN("could not find the samtools source sub-directory; some commands disabled")
fi

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
