#!/usr/bin/env python
# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved

import os
import sys
from optparse import OptionParser
import progressbar

fields = ['qname', 'flag', 'rname', 'pos', 'mapq',
          'cigar', 'rnext', 'pnext', 'tlen', 'seq',
          'qual', 'rg', 'pg', 'md', 'nm', 'as', 'fz', 
          'xa', 'xs', 'xt']

class Record(object):
    unNamedRange = (0, 10)
    def __init__(self, tokens):
        self._set_attrs( tokens )
        """
        self.fields_len = len(tokens)
        for x in xrange(len(fields)):
            try:
                setattr(self, fields[x], tokens[x])
            except:
                print tokens
                print fields
                print "len(tokens)=", len(tokens)
                print "len(fields)=", len(fields)
                sys.exit(1)
        """
    def __str__(self):
        return '\t'.join( [ getattr(self,fields[x]) for x in xrange(self.fields_len) ] )

    def _set_attrs( self, tokens ):
        for x in xrange( self.unNamedRange[1] ):
            setattr(self, fields[x], tokens[x])

        for x in xrange( self.unNamedRange[1] + 1, len( tokens ) ):
            fieldName = tokens[x].split(':')[0].lower()
            setattr(self, fieldName, tokens[x])


class Sam(object):

    def __init__(self, sam, full_qname ):
        self.records = {}
        self.sam = sam
        self._parse_file( full_qname )

    def _parse_file(self, full_qname ):
        fp = open( self.sam, "r" )
        fp = fp.readlines()
        widgets = [
                    "parsing %s: " % ( self.sam ),
                    progressbar.Percentage(),
                    progressbar.Bar()
                    ]
        pbar = progressbar.ProgressBar(widgets=widgets, maxval=len(fp)).start()
        for i,line in enumerate(fp):
            if line[0] == '@':
                #not supporting headers for now
                continue
            line = line.rstrip()
            tokens = line.split('\t')
            rec = Record(tokens)
            self.records[ self._hash_name( rec.qname, full_qname ) ] = rec
            pbar.update(i+1)
        pbar.finish()

    def _hash_name( self, name, full_qname ):
        if full_qname:
            return name
        else:
            return ':'.join(name.split(':')[1:])

def diff_field(field1, field2):
    """returns true if field1 == field2"""
    return field1 == field2

def main(options):
    sam1 = Sam(options.sam1, options.full_qname)
    sam2 = Sam(options.sam2, options.full_qname)
    fields = options.fields.split(',')
    widgets = [
                "diffing: ", 
                progressbar.Percentage(),
                progressbar.Bar()
                ]
    pbar = progressbar.ProgressBar( widgets=widgets, maxval=len( sam1.records) ).start()

    for i, read in enumerate( (sam1.records.keys() ) ):
        if read not in sam2.records:
            print "read", read,"not in sam2"
        else:
            if 0 < options.min_mapq:
                mapq1 = geattr( sam1.records[ read ], 'mapq' )
                mapq2 = geattr( sam2.records[ read ], 'mapq' )
                if mapq1 < options.min_mapq and mapq2 < options.min_mapq:
                    continue
            diff_str = "[%s]" % (read)
            for field in fields:
                cont = [False, False] #var to track these 2 exceptions below to see which one failed
                try:
                    attr1 = getattr( sam1.records[ read ], field )
                except:
                    #print sam1.sam, read, "doesn't have the: ", field, " tag"
                    cont[0] = True
                try:
                    attr2 = getattr( sam2.records[ read ], field )
                except:
                    cont[1] = True
                if cont[0] and cont[1]:
                    continue
                elif cont[0] and not cont[1]:
                    print read, sam1.sam, "has the field: ", field, " and", sam2.sam, "does not"
                    continue
                elif not cont[0] and cont[1]:
                    print read, sam2.sam, "has the field: ", field, " and", sam1.sam, "does not"
                    continue
                if not diff_field( getattr(sam1.records[ read ], field), getattr(sam2.records[ read ], field) ):
                    diff_str = "%s -- %s[%s]=%s %s[%s]=%s" % (diff_str, sam1.sam, field, str(attr1), sam2.sam, field, str(attr2))
            if len(diff_str) > len(read) + 2:
                print diff_str
        pbar.update(i+1)
    pbar.finish()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--sam1', help="first sam file to diff.  will be called sam1 in diff out", dest='sam1')
    parser.add_option('--sam2', help="second sam file to diff.  will be called sam2 in diff out", dest='sam2')
    parser.add_option('--fields',
                      help="comma seperated list of fields:%s\t\t\t\t\t to diff between the"
                      "sam records use names from the same spec. for optional tags"
                      "use their 2 letter abbreviation.  Default:  pos" % (str(fields)),
                      dest='fields', default=['pos'])
    parser.add_option('--full-qname', help="keep the full query name", dest='full_qname', action="store_true", default=False)
    parser.add_option('--min-mapq', help="examine only those records with a given minimum mapping quality", dest="min_mapq", default=0)
    options, args = parser.parse_args()
    main(options)
