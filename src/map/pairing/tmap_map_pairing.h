/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_MAP_PAIRING_H
#define TMAP_MAP_PAIRING_H

// TODO
enum {
    TMAP_MAP_PAIRING_SAME_STRAND=0,
    TMAP_MAP_PAIRING_OPPOSITE_STRAND=1
};

// TODO
enum {
    TMAP_MAP_PAIRING_POSITIONING_AB=0,
    TMAP_MAP_PAIRING_POSITIONING_BA,
    TMAP_MAP_PAIRING_POSITIONING_NONE
};

// TODO
void
tmap_map_pairing_pick_pairs(tmap_map_sams_t *one, tmap_map_sams_t *two, tmap_seq_t *one_seq, tmap_seq_t *two_seq,
                            tmap_rand_t *rand, tmap_map_opt_t *opt);

#endif
