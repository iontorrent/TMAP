#!/usr/bin/env python
# Copyright (C) 2011 Ion Torrent Systems, Inc. All Rights Reserved

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
        qseq = qseq[:0] + '' + qseq[0+count:]
    elif operation == 4 or operation == 1:
        qseq = qseq[:0] + '' + qseq[0+count:]

  start = 0
  result = ""
  md_split = re.findall('(\d+)|\^([gatcnryswkmbdhv]+)|([gatcnryswkmbdhvn]+)', md, re.IGNORECASE)
  for i in md_split:
    if i[0]:
      end = start + int(i[0])
      result += seq[start:end]
      start += int(i[0])
    elif i[1]:
      result += i[1]
    elif i[2]:
      result += i[2]
      start += len(i[2])
  return result

def translate(ref, read):
    read = read[::-1]
    read = read.translate(string.maketrans('ACGTacgt', 'TGCAtgca'))
    ref = ref[::-1]
    ref = ref.translate(string.maketrans('ACGTacgtRYKMBVDH', 'TGCAtgcaYRMKVBHD'))
    return (ref, read)

def left_justify(read, ref, matcha):

  #print "before"
  #print read
  #print matcha
  #print ref
  read = list(read)
  ref = list(ref)
  matcha = list(matcha)

  prev_del = 0;
  prev_ins = 0;
  start_ins = 0;
  start_del = 0;
  end_ins = 0;
  end_del = 0;
  justified = 0;

  for i in range(len(matcha)):
    if read[i] == '-':
      if prev_del == 0:
        start_del = i
      prev_del = 1
      end_del = i
      prev_ins = 0
      start_ins = end_ins = -1
      i += 1

    elif ref[i] == '-':
      if prev_ins == 0:
        start_ins = i
      prev_ins = 1
      end_ins = i
      prev_del = 0
      start_del = -1
      end_del = -1
      i += 1

    else:
      if prev_del == 1:
        #print "here"
        start_del -= 1
        while start_del >= 0 and read[start_del] != '-' and ref[start_del] != '-' and ref[start_del] == ref[end_del]:
          #print "inside"
          c = read[end_del]
          read[end_del] = read[start_del]
          read[start_del] = c
          c = matcha[end_del]
          matcha[end_del] = matcha[start_del]
          matcha[start_del] = c
          start_del -= 1
          end_del -= 1
          justified = 1
        end_del += 1
        i = end_del 

      elif prev_ins == 1:
        start_ins -= 1
        while start_ins >= 0 and read[start_ins] != '-' and ref[start_ins] != '-' and read[start_ins] == read[end_ins]:
          c = ref[end_ins]
          ref[end_ins] = ref[start_ins]
          ref[start_ins] = c
          c = matcha[end_ins]
          matcha[end_ins] = matcha[start_ins]
          matcha[start_ins] = c
          start_ins -= 1
          end_ins -= 1
          justified = 1
        end_ins += 1
        i = end_ins

      else:
        pass

      prev_del = prev_ins = 0;
      start_del = start_ins = end_del = end_ins = -1
 
  read = ''.join(read)
  ref = ''.join(ref)
  matcha = ''.join(matcha)
  return read, ref, matcha
      
def alignment(strand, md, cigar, qseq):
  ref, read, match = ("", "", "")
  qref_i, qseq_i, n = (0, 0, 0)
  iupac = {'R' : 'AG', 'Y' : 'CT', 'S' : 'GC', 'W' : 'AT', 'K' : 'GT', 'M' : 'AC', 'B' : 'CGT', 'D' : 'AGT', 'H' : 'ACT', 'V' : 'ACG', 'N' : 'ACGT',}

  qref = get_qref(md, cigar, qseq)
  qref = qref.upper()
  qseq = qseq.upper()

  if strand == 1:
    (qref, qseq) = translate(qref, qseq)
    cigar = reversed(cigar)
  #make the alignment
  for i in cigar:
    #print i
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
    elif operation == 1:
      read += qseq[qseq_i:qseq_i+count]
      qseq_i += count
      ref += "-" * count
      match += "+" * count
    elif operation == 2:
      read += "-" * count
      match += "-" * count
      ref += qref[qref_i:qref_i+count]
      qref_i += count
    elif operation == 4:
      #read += qseq[0:count]
      #match += "S" * count
      #ref += "-" * count
      qseq = qseq[count:]
        #print "readh = " + read
            #print "match = " + match
            #print "refdh = " + ref
    n+=1
    #if strand == 1:
    #read, ref, match = left_justify(read, ref, match)

  #print read
  #print ref
  #print match
  #if strand == 1:
  #  read1, ref1, match1 = left_justify(read, ref, match)
  #  return (read1.upper(), ref1.upper(), match1)
  #else:
  return (read.upper(), ref.upper(), match)

def main(options):
  if re.search(r"sam$", options.sam):
    sam = pysam.Samfile(options.sam,"r")
  elif re.search(r"bam$", options.sam):
    sam = pysam.Samfile(options.sam,"rb")

  count = 0
  for lines in sam.fetch():
    if lines.is_unmapped:
      continue
    count+=1
    if (count%10000) == 0:
        sys.stderr.write(str(count)+"\r")
    
    qname = lines.qname
    md = lines.opt('MD')
    cigar = lines.cigar
    qseq = lines.seq
    strand = 0
    if lines.is_reverse:
      strand = 1
    
    qdna, tdna, matcha = alignment(strand, md, cigar, qseq)
    #print qdna
    #print matcha
    #print tdna

if __name__ == '__main__':
  parser = OptionParser()
  parser.add_option('--sam', help="input sam file", dest='sam')
  if len(sys.argv[1:]) < 1:
    parser.print_help()
  else:
    options, args = parser.parse_args()
    main(options)
