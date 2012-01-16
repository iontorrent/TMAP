/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_TIME_H
#define TMAP_TIME_H

/*! 
  CPU and Realtime timing
  */

/*!
  @return returns the elapsed CPU time.
 */
double 
tmap_time_cputime();

/*!
  @return returns the real time.
 */
double 
tmap_time_realtime();

#endif
