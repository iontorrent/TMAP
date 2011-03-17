/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_ALL_H
#define TMAP_MAP_ALL_H

#include <sys/types.h>

/*!
  Parses the command line options and stores them in the options structure
  @param  argc  the number of arguments
  @param  argv  the argument list
  @param  opt   pointer to the options
  @return       1 if successful, 0 otherwise
  */
int32_t
tmap_map_all_opt_parse(int argc, char *argv[], tmap_map_opt_t *opt);

/*! 
  main-like function for 'tmap map_all'
  @param  argc  the number of arguments
  @param  argv  the argument list
  @return       0 if executed successful
  */
int 
tmap_map_all_main(int argc, char *argv[]);

#endif 
