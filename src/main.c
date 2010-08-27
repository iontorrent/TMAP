#include <stdio.h>
#include <string.h>
#include <config.h>
#include <stdint.h>

#include "fmap_error.h"

extern int fmap_index(int argc, char *argv[]);
extern int fmap_exact(int argc, char *argv[]);

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
	fprintf(stderr, "         index\n");
        fprintf(stderr, "Debugging:\n");
	fprintf(stderr, "         exact\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if(argc < 2) return usage();
	else if (0 == strcmp("index", argv[1])) return fmap_index(argc-1, argv+1);
	else if (0 == strcmp("exact", argv[1])) return fmap_exact(argc-1, argv+1);
	else {
            fmap_error1(PACKAGE, argv[1], Exit, CommandLineArgument);
	}
	return 0;
}

