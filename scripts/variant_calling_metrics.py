#!/usr/bin/env python
# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved

import os
import re
import sys
from optparse import OptionParser
from itertools import *
import progressbar

def diff_int(val1, val2, wiggle):
    """returns true if abs(val2 - val1) <= wiggle"""
    if val1 < val2:
        return val2 - val1 <= wiggle
    else:
        return val1 - val2 <= wiggle

def diff_field(field1, field2):
    """returns true if field1 == field2"""
    return field1 == field2

def aln_len(cigar1, cigar2):
    """returns alignment length"""
    aln1 = sum(int(x) for x in cigar1.split("|")[0:-1])
    aln2 = sum(int(x) for x in cigar2.split("|")[0:-1])
    return aln1, aln2

def main(options):
    sam1 = options.sam1
    sam2 = options.sam2

    reads_mapped_with_min_mapq_only_in_sam1 = 0
    reads_mapped_with_min_mapq_only_in_sam2 = 0
    reads_mapped_with_min_mapq_in_sam1_and_less_min_mapq_in_sam2 = 0
    reads_mapped_only_in_sam1 = 0
    reads_mapped_only_in_sam2 = 0
    reads_mapped_with_min_mapq_diff_aln_pos = 0
    reads_mapped_diff_aln_pos = 0
    reads_mapped_with_min_mapq_le_aln_len = 0
    reads_mapped_with_min_mapq_ge_aln_len = 0
    reads_mapped_le_aln_len = 0
    reads_mapped_ge_aln_len = 0
    reads_mapped_with_min_mapq_diff_cigar_string = 0
    reads_mapped_diff_cigar_string = 0

    count = 0
    for line1, line2 in izip(open(sam1, "r"), open(sam2, "r")):
        count+=1
        if (count%10000) == 0:
            sys.stderr.write(str(count)+"\r")

        if line1.startswith("@") and line2.startswith("@"):
            continue
        #print line1, line2
        temp_line1 = line1.split()
        temp_line2 = line2.split()

        if temp_line1[0].strip() == temp_line2[0].strip():
            flag_sam1 = int(temp_line1[1])
            mapq_sam1 = int(temp_line1[4])
            rname_sam1 = temp_line1[2]
            pos_sam1 = int(temp_line1[3])
            cigar_sam1 = temp_line1[5]

            flag_sam2 = int(temp_line2[1])
            mapq_sam2 = int(temp_line2[4])
            rname_sam2 = temp_line2[2]
            pos_sam2 = int(temp_line2[3])
            cigar_sam2 = temp_line2[5]

            if (flag_sam1 & 4) == 0 and (flag_sam2 & 4) == 0: #both reads mapped
                if not diff_field( rname_sam1, rname_sam2 ): #if chr are not equal
                    reads_mapped_diff_aln_pos += 1
                elif not diff_int( pos_sam1, pos_sam2, options.wiggle ): #if pos are not same
                    reads_mapped_diff_aln_pos += 1
                else:
                    (cigar_sam1_MD, count1) = re.subn(r'[MD]', '|', re.subn(r'\d+[SIHNP]', '', cigar_sam1)[0])
                    (cigar_sam2_MD, count2) = re.subn(r'[MD]', '|', re.subn(r'\d+[SIHNP]', '', cigar_sam2)[0])
                    (sam1_aln_len, sam2_aln_len) = aln_len(cigar_sam1_MD, cigar_sam2_MD)
                    if sam1_aln_len < sam2_aln_len:
                        reads_mapped_le_aln_len += 1
                    elif sam1_aln_len > sam2_aln_len:
                        reads_mapped_ge_aln_len += 1

                    if not diff_field (cigar_sam1, cigar_sam2):
                        reads_mapped_diff_cigar_string += 1

                if mapq_sam1 > options.min_mapq and mapq_sam2 > options.min_mapq: #if mapq > options.mapq
                    if not diff_field( rname_sam1, rname_sam2 ): #if chr are not equal
                        reads_mapped_with_min_mapq_diff_aln_pos += 1
                    elif not diff_int( pos_sam1, pos_sam2, options.wiggle ): #if pos are not same
                        reads_mapped_with_min_mapq_diff_aln_pos += 1
                    else:
                        (cigar_sam1_MD, count1) = re.subn(r'[MD]', '|', re.subn(r'\d+[SIHNP]', '', cigar_sam1)[0])
                        (cigar_sam2_MD, count2) = re.subn(r'[MD]', '|', re.subn(r'\d+[SIHNP]', '', cigar_sam2)[0])
                        (sam1_aln_len, sam2_aln_len) = aln_len(cigar_sam1_MD, cigar_sam2_MD)
                        if sam1_aln_len < sam2_aln_len:
                            reads_mapped_with_min_mapq_le_aln_len += 1
                        elif sam1_aln_len > sam2_aln_len:
                            reads_mapped_with_min_mapq_ge_aln_len += 1
                            reads_mapped_ge_aln_len += 1

                        if not diff_field (cigar_sam1, cigar_sam2):
                            reads_mapped_with_min_mapq_diff_cigar_string += 1

            elif (flag_sam1 & 4) == 0 and (flag_sam2 & 4) != 0: #read2 is unmapped 
                reads_mapped_only_in_sam1 += 1
                if mapq_sam1 > options.min_mapq:
                    reads_mapped_with_min_mapq_only_in_sam1 += 1
            elif (flag_sam1 & 4) != 0 and (flag_sam2 & 4) == 0: #read1 is unmapped
                reads_mapped_only_in_sam2 += 1
                if mapq_sam2 > options.min_mapq:
                    reads_mapped_with_min_mapq_only_in_sam2 += 1

        else:
            print "Sam record " + temp_line1[0] + " from sam file " + sam1 + " not the same as " + temp_line2[0] + " from sam file " + sam2

    #print options.ver1 + "_" + options.ver2 + "_" + options.sam1 + "_" + options.sam2
    print "Reads mapping only in version1 = " + str(reads_mapped_only_in_sam1)
    print "Reads mapping only in version2 = " + str(reads_mapped_only_in_sam2)
    print "Reads above mapqv " + str(options.min_mapq) + " mapping only in version1 = " + str(reads_mapped_with_min_mapq_only_in_sam1)
    print "Reads above mapqv " + str(options.min_mapq) + " mapping only in version2 = " + str(reads_mapped_with_min_mapq_only_in_sam2)
    print ""
    print "Reads mapping in both versions:"
    print "At different alignment positions = " + str(reads_mapped_diff_aln_pos)
    print "With different alignment lengths = " + str(reads_mapped_le_aln_len + reads_mapped_ge_aln_len)
    print "With version1 smaller alignment lengths = " + str(reads_mapped_le_aln_len)
    print "With version2 smaller alignment lengths = " + str(reads_mapped_ge_aln_len)
    print "With different cigar strings = " + str(reads_mapped_diff_cigar_string)
    print ""
    print "Reads mapping in both versions above mapq " + str(options.min_mapq) + " :" 
    print "At different alignment positions above mapq " + str(options.min_mapq) + " = " + str(reads_mapped_with_min_mapq_diff_aln_pos)
    print "With different alignment lengths above mapq " + str(options.min_mapq) + " = " + str(reads_mapped_with_min_mapq_le_aln_len + reads_mapped_with_min_mapq_ge_aln_len)
    print "With version1 smaller alignment lengths above mapq " + str(options.min_mapq) + " = " + str(reads_mapped_with_min_mapq_le_aln_len)
    print "With version2 smaller alignment lengths above mapq " + str(options.min_mapq) + " = " + str(reads_mapped_with_min_mapq_ge_aln_len)
    print "With different cigar strings above mapq " + str(options.min_mapq) + " = " + str(reads_mapped_with_min_mapq_diff_cigar_string)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--sam1', help="first sam file to diff.  will be called sam1 in diff out", dest='sam1')
    parser.add_option('--sam2', help="second sam file to diff.  will be called sam2 in diff out", dest='sam2')
    parser.add_option('--wiggle', help="allow positions to be within this window to be considered the same", default=0, dest='wiggle')
    #parser.add_option('--ver1', help="tmap version 1", dest='ver1')
    #parser.add_option('--ver2', help="tmap version 2", dest='ver2')
    parser.add_option('--min-mapq', help="the minimum mapping quality for secondary analysis", type="int", dest="min_mapq", default=8)
    if len(sys.argv[1:]) < 1:
        parser.print_help()
    else:
        options, args = parser.parse_args()
        main(options)
