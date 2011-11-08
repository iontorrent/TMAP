import os
import sys
from optparse import OptionParser

fields = ['qname', 'flag', 'rname', 'pos', 'mapq',
          'cigar', 'rnext', 'pnext', 'tlen', 'seq',
          'qual', 'rg', 'md', 'nm', 'as', 'xs', 'xt']

class Record(object):
    
    def __init__(self, tokens):
        for x in xrange(len(fields)):
            setattr(self, fields[x], tokens[x])

class Sam(object):
    
    def __init__(self, sam):
        self.records = {}
        self.sam = sam
        
    def _parse_file(self):
        fp = open(sam, "r")
        for line in fp:
            if line[0] == '@':
                #not supporting headers for now
                continue
            line = line.rstrip()
            tokens = line.split('\t')
            rec = Record(tokens)
            self.records[ rec.qname ] = rec

            
            
def main(options):
    sam1 = Sam(options.sam1)
    sam2 = Sam(options.sam2)
    
    for read in (sam1.records.keys()):
        if read not in sam2.records:
            print "read ", read, "not in sam2"
            
    for read in (sam2.records.keys()):
        if read not in sam1.records:
            print "read ", read, "not in sam1"
            
    

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--sam1', help="first sam file to diff.  will be called sam1 in diff out", dest='sam1')
    parser.add_option('--sam2', help="second sam file to diff.  will be called sam2 in diff out", dest='sam2')    
    options, args = parser.parse_args()
    
    main(options)
    