#!/usr/bin/env python
# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved

import sys
from optparse import OptionParser

class Fastq:
    name = ''
    seq = ''
    comment = ''
    qual = ''
    fh = None
    eof = False
    loaded = False

    def __init__(self, fn):
        self.fh = open(fn, 'r')
        self.loaded = False
        self.eof = False

    def read(self):
        self.name = self.fh.readline().rstrip()
        if '' == self.name:
            self.eof = True
            self.loaded = False
            self.fh.close()
        else:
            self.seq = self.fh.readline().rstrip()
            self.comment = self.fh.readline().rstrip()
            self.qual = self.fh.readline().rstrip()
            s = self.name.split(':')
            self.row = int(s[1]) 
            self.col = int(s[2])
            self.eof = False
            self.loaded = True

    def write(self, fh):
        fh.write(self.name + '\n' + self.seq + '\n' + self.comment + '\n' + self.qual + '\n')
        self.loaded = False

def main(options):
    fh_in_two = open(options.fn_read_two, 'r')

    f_one = Fastq(options.fn_read_one)
    f_two = Fastq(options.fn_read_two)

    fh_out_one = open(options.fn_prefix + '.read1.fastq', 'w')
    fh_out_two = open(options.fn_prefix + '.read2.fastq', 'w')
    fh_out_single = open(options.fn_prefix + '.single.fastq', 'w')


    while (not f_one.eof) or (not f_two.eof) :
        if (not f_one.loaded) and (not f_one.eof):
            f_one.read()
        if (not f_two.loaded) and (not f_two.eof):
            f_two.read()
        
        case = -1
        if (not f_one.loaded) and (not f_two.loaded):
            case = -1
        elif not f_one.loaded:
            case = 2
        elif not f_two.loaded:
            case = 0
        elif (f_one.row < f_two.row) or (f_one.row == f_two.row and f_one.col < f_two.col):
            case = 0
        elif f_one.row == f_two.row and f_one.col == f_two.col:
            case = 1
        else:
            case = 2

        if 0 == case:
            f_one.write(fh_out_single)
        elif 1 == case:
            f_one.write(fh_out_one)
            f_two.write(fh_out_two)
        elif 2 == case:
            f_two.write(fh_out_single);

    fh_out_one.close()
    fh_out_two.close()
    fh_out_single.close()

def check_option(parser, value, name):
    if None == value:
        print 'Option ' + name + ' required.\n'
        parser.print_help()
        sys.exit(1)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-1', '--fn-read-one', dest='fn_read_one', default=None, help="Read one")
    parser.add_option('-2', '--fn-read-two', dest='fn_read_two', default=None, help="Read two")
    parser.add_option('-p', '--fn-prefix', dest='fn_prefix', default=None, help="Output prefix")

    (options, args) = parser.parse_args()
    if len(args) != 0:
        parser.print_help()
        sys.exit(1)

    check_option(parser, options.fn_read_one, '-1')
    check_option(parser, options.fn_read_two, '-2')
    check_option(parser, options.fn_prefix, '-p')

    main(options)
