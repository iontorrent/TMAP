import os
import sys
from optparse import OptionParser

fields = ['qname', 'flag', 'rname', 'pos', 'mapq',
          'cigar', 'rnext', 'pnext', 'tlen', 'seq',
          'qual', 'rg', 'pg', 'md', 'nm', 'as', 'nh', 
          'xs', 'xt']

class Record(object):
    
    def __init__(self, tokens):
        self.fields_len = len(tokens)
        for x in xrange(len(tokens)):
            setattr(self, fields[x], tokens[x])
    def __str__(self):
        return '\t'.join( [ getattr(self,fields[x]) for x in xrange(self.fields_len) ] )

class Sam(object):
    
    def __init__(self, sam):
        self.records = {}
        self.sam = sam
        self._parse_file()
        
    def _parse_file(self):
        fp = open( self.sam, "r" )
        for line in fp:
            if line[0] == '@':
                #not supporting headers for now
                continue
            line = line.rstrip()
            tokens = line.split('\t')
            rec = Record(tokens)
            self.records[ rec.qname ] = rec

            

def diff_field(field1, field2):
    """returns true if field1 == field2"""
    return field1 == field2

def main(options):
    sam1 = Sam(options.sam1)
    sam2 = Sam(options.sam2)
    fields = options.fields.split(',')
    for read in (sam1.records.keys()):
        if read not in sam2.records:
            print "read", read,"not in sam2"
        else:
            for field in fields:
                
                try:
                    attr1 = getattr( sam1.records[ read ], field )
                except:
                    print sam1.sam, read, "doesn't have the: ", field, " tag"
                    continue
                try:
                    attr2 = getattr( sam2.records[ read ], field )
                except:
                    print sam2.sam, read, "doesn't have the: ", field, " tag"
                    continue
                if not diff_field( getattr(sam1.records[ read ], field), getattr(sam2.records[ read ], field) ):
                    print "Different at field: ", field
                    print sam1.sam, getattr(sam1.records[ read ], field)
                    print sam2.sam, getattr(sam2.records[ read ], field)
                    print sam1.sam,": ", str(sam1.records[read])
                    print sam2.sam,": ", str(sam2.records[read])

    

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--sam1', help="first sam file to diff.  will be called sam1 in diff out", dest='sam1')
    parser.add_option('--sam2', help="second sam file to diff.  will be called sam2 in diff out", dest='sam2')
    parser.add_option('--fields', help="""comma seperated list of fields to diff between the sam records.  
                                          use names from the same spec.  for optional tags, use their 2 letter 
                                          abbreviation""", dest='fields', default=[])
    options, args = parser.parse_args()
    
    main(options)
    
