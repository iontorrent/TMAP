#!/usr/bin/env python
# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved

import os
import re
import sys
from optparse import OptionParser
import pysam
import string

#MATCH     = 0 
#INS       = 1 
#DEL       = 2 
#REF_SKIP  = 3 
#SOFT_CLIP = 4 
#HARD_CLIP = 5 
#PAD       = 6

def get_qref(md, cigar, qseq):
  seq = ""
  for i in cigar:
    operation, count = i[0], int(i[1])
    #print operation, count
    if operation == 0:
        seq += qseq[0:count]
        #qseq = qseq.replace(qseq[0:count], '')
        qseq = qseq[:0] + '' + qseq[0+count:]
    elif operation == 4 or operation == 1:
        #print "beforeeeeeeeeeeeeeee = " + qseq
        #qseq = qseq.replace(qseq[0:count], '')
        qseq = qseq[:0] + '' + qseq[0+count:]
        #print "qseqqqqqqqqqqqqqqqqqqqqqqqqqqqq = " + qseq

    #print "seq = " + seq
  start = 0
  result = ""
  md_split = re.findall('(\d+)|\^([gatcnryswkmbdhv]+)|([gatcnryswkmbdhvn]+)', md, re.IGNORECASE)
  #print md_split
  for i in md_split:
    #print i
    if i[0]:
      end = start + int(i[0])
      result += seq[start:end]
      start += int(i[0])
      #print result
    elif i[1]:
      result += i[1]
      #print result
    elif i[2]:
      result += i[2]
      start += len(i[2])
      #print result
    #print "result = " + result 
  return result

def translate(read, ref, match):
    read = read[::-1]
    read = read.translate(string.maketrans('ACGTacgt', 'TGCAtgca'))
    ref = ref[::-1]
    ref = ref.translate(string.maketrans('ACGTacgtRYKMBVDH', 'TGCAtgcaYRMKVBHD'))
    match = match[::-1]
    return (read, ref, match)

def alignment(strand, md, cigar, qseq):
  ref, read, match = ("", "", "")
  qref_i, qseq_i, n = (0, 0, 0)
  iupac = {'R' : 'AG', 'Y' : 'CT', 'S' : 'GC', 'W' : 'AT', 'K' : 'GT', 'M' : 'AC', 'B' : 'CGT', 'D' : 'AGT', 'H' : 'ACT', 'V' : 'ACG', 'N' : 'ACGT',}

  qref = get_qref(md, cigar, qseq)
  #print qref
  qref = qref.upper()
  qseq = qseq.upper()

  #make the alignment
  #print "cigar = " + str(cigar)
  for i in cigar:
    #print i[0]
    operation, count = i[0], int(i[1])
    if operation == 0:
      for i in range(0, count):
        if qref[qref_i+i : qref_i+i+1] in iupac:
          if qseq[qseq_i+i : qseq_i+i+1] in iupac[qref[qref_i+i : qref_i+i+1]]:
            match += "|"
          else:
            match += " "
        elif qseq[qseq_i+i:qseq_i+i+1] == qref[qref_i+i:qref_i+i+1]:
          match += "|"
        else:
          match += " "
      read += qseq[qseq_i:qseq_i+count]
      qseq_i += count
      ref += qref[qref_i:qref_i+count]
      qref_i += count
        #print "readh = " + read
        #print "match = " + match
        #print "refdh = " + ref
    elif operation == 1:
      read += qseq[qseq_i:qseq_i+count]
      qseq_i += count
      ref += "-" * count
      match += "-" * count
        #print "readh = " + read
        #print "match = " + match
        #print "refdh = " + ref
    elif operation == 2:
      read += "-" * count
      match += "-" * count
      ref += qref[qref_i:qref_i+count]
      qref_i += count
              #print "readh = " + read
            #print "match = " + match
            #print "refdh = " + ref
    elif operation == 4:
      read += qseq[0:count]
      match += "S" * count
      ref += "-" * count
      qseq = qseq[count:]
        #print "readh = " + read
            #print "match = " + match
            #print "refdh = " + ref
    n+=1

  if strand == 1:
    #print "readh = " + read
    #print "match = " + match
    #print "refdh = " + ref
    (read, ref, match) = translate(read, ref, match)

  return (read.upper(), ref.upper(), match)


class HPStats:

    def __init__(self, maxhp, maxrl):
        self._data = list()
        self._maxhp = maxhp
        self._maxrl = maxrl

    def add(self, l, i, aln_type):
        ''' l is the hp length (one-based), i is the read index (zero-based), aln_type is (0 - non-indel, 1 - insertion, 2 - deletion)  '''
        if self._maxhp < l:
            l = self._maxhp
        if self._maxrl <= i:
            return # do not add
        while len(self._data) < l:
            self._data.append(list())
        while len(self._data[l-1]) <= i:
            self._data[l-1].append([0,0,0])
        self._data[l-1][i][0] += 1 # total
        if 0 < aln_type:
            self._data[l-1][i][aln_type] += 1

    def get(self, l, i, aln_type):
        if self._maxhp < l:
            l = self._maxhp
        if len(self._data) < l:
            return None
        if len(self._data[l-1]) <= i:
            return None
        if aln_type < 0 or 2 < aln_type:
            return None
        return self._data[l-1][i][aln_type]

    def dump(self):
        for i in xrange(len(self._data)): # hp length
            for j in xrange(len(self._data[i])): # read index
                print "%d\t%d\t%d\t%d\t%d" % (i+1, j+1, \
                        self._data[i][j][0], \
                        self._data[i][j][1], \
                        self._data[i][j][2])

    def dump_table(self, num_reads):
        num_total = num_insertions = num_deletions = 0
        for i in xrange(len(self._data)): # hp length
            for j in xrange(len(self._data[i])): # read index
                num_total += self._data[i][j][0]
                num_insertions += self._data[i][j][1]
                num_deletions += self._data[i][j][2]
        print "%d\t%d\t%.2f\t%.2f" % (num_insertions, num_deletions, \
                100.0 * (num_insertions + num_deletions) / num_total, \
                (num_insertions + num_deletions) / num_reads)


def main(options):
  if re.search(r"sam$", options.sam):
    sam = pysam.Samfile(options.sam,"r")
  elif re.search(r"bam$", options.sam):
    sam = pysam.Samfile(options.sam,"rb")

  if options.hpstats:
      hpstats = HPStats(options.hpstats_maxhp, options.hpstats_maxrl)

  n = 0
  for read in sam.fetch():
    if read.is_unmapped:
        continue
    n += 1
    md = read.opt('MD')
    cigar = read.cigar
    qseq = read.seq
    strand = 0
    if read.is_reverse:
      strand = 1

    #print str(strand), md, cigar, qseq
    read, ref, match = alignment(strand, md, cigar, qseq)    
    if not options.hpstats:
        print read
        print match
        print ref
    else:

        read_i = [0 for i in xrange(len(match))]
        j = 0
        for i in xrange(len(match)):
            read_i[i] = j
            if '-' != read[i]:
                j += 1

        # init
        if 1 == strand:
            (read, ref, match) = translate(read, ref, match)
                
        i = 0
        while i < len(match) and 'S' == match[i]:
            i += 1

        i_end = i
        while i_end < len(match) and 'S' != match[i_end]:
            i_end += 1

        while i < i_end:
            # go through an hp
            start = end = i

            read_l = 0
            ref_l = 0

            # find the end
            if '-' == read[end]:
                end = start;
                while end < i_end and ('-' == read[end] or ref[start] == read[end]) and ref[start] == ref[end]:
                    ref_l += 1
                    if '-' != read[end]:
                        read_l += 1
                    end += 1
            else:
                end = start
                while end < i_end and ('-' == read[end] or read[start] == read[end]) and ('-' == ref[end] or read[start] == ref[end]):
                    if '-' != read[end]:
                        read_l += 1
                    if '-' != ref[end]:
                        ref_l += 1
                    end += 1

            if 0 < read_l or 0 < ref_l:
                if options.hpstats_verbose:
                    print "HP Found at:%d ref:%d read:%d" % (read_i[start], ref_l, read_l) 
                if read_l < ref_l: # deletion
                    hpstats.add(ref_l, read_i[start], 2)
                elif ref_l < read_l: # insertion
                    if options.hpstats_ins_by_ref:
                        hpstats.add(ref_l, read_i[start], 1)
                    else:
                        hpstats.add(read_l, read_i[start], 1)
                else:
                    # normal
                    hpstats.add(ref_l, read_i[start], 0)

            if end == start:
                i += 1
            else:
                i = end
        if options.hpstats_verbose:
            print read
            print match
            print ref
  if options.hpstats:
      if options.hpstats_table:
          hpstats.dump_table(n)
      else:
          hpstats.dump()

if __name__ == '__main__':
  parser = OptionParser()
  parser.add_option('--sam', help="input sam file", dest='sam', default=None)
  parser.add_option('--hpstats', help="hpstats", action="store_true", dest='hpstats', default=False)
  parser.add_option('--hpstats-ins-by-ref', help="use the reference hp length for insertions", action="store_true", dest='hpstats_ins_by_ref', default=False)
  parser.add_option('--hpstats-maxhp', type=int, help="maximum homopolymer length for hpstats", dest='hpstats_maxhp', default=9)
  parser.add_option('--hpstats-maxrl', type=int, help="maximum read length for hpstats", dest='hpstats_maxrl', default=100)
  parser.add_option('--hpstats-table', help="dump indel summary for hpstats", dest='hpstats_table', action="store_true", default=False)
  parser.add_option('--hpstats-verbose', help="hpstats verbose", dest='hpstats_verbose', action="store_true", default=False)
  if len(sys.argv[1:]) < 1:
    parser.print_help()
  else:
    options, args = parser.parse_args()
    main(options)
