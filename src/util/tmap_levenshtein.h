/* Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved */
#ifndef TMAP_LEVENSHTEIN_H
#define TMAP_LEVENSHTEIN_H

/*!
 This function implements the Damerau-Levenshtein algorithm to calculate a distance between strings.
 @param  string1                string 1
 @param  string2                string 2
 @param  swap_penalty           the swap penalty
 @param  substitution_penalty   the substitution penalty
 @param  insertion_penalty      the insertion penalty
 @param  deletion_penalty       the deletion penalty
 @return  the edit distance
 */
int 
tmap_levenshtein(const char *string1, const char *string2,
	int swap_penalty, int substitution_penalty,
	int insertion_penalty, int deletion_penalty);

#endif
