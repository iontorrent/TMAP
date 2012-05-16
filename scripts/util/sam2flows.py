#!/usr/bin/env python

import optparse
import random
import sys

def reverse_compliment(seq):
    c = {'a':'t', 'c':'g', 'g':'c', 't':'a', 'n':'n',
            'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    return "".join([c.get(nt, '') for nt in seq[::-1]])

def make_flow_signal(signal):
    r = random.gauss(0, 0.25)
    if 0 == signal:
        while r < 0 or 0.5 < r:
            r = random.gauss(0, 0.25)
    else: 
        while r < -0.5 or 0.5 < r:
            r = random.gauss(0, 0.25)
    return (100 * signal) + int(r * 100)

def check_option(parser, value, name):
    if None == value:
        print 'Option ' + name + ' required.\n'
        parser.print_help()
        sys.exit(1)

def __main__():
    parser = optparse.OptionParser()

    parser.add_option( '-f', '--flow-order', dest='flow_order', help='the flow order to simulate' )
    parser.add_option( '-k', '--key-sequence', dest='key_sequence', help='the key sequence to simulate' )
    parser.add_option( '-s', '--sam-file', dest='sam_file', help='the input SAM file' )
    parser.add_option( '-S', '--seed', dest='seed', help='the seed for the random number generator' )

    (options, args) = parser.parse_args()
    
    if len(args) != 0:
        parser.print_help()
        sys.exit(1)
    check_option(parser, options.flow_order, '-f')
    check_option(parser, options.key_sequence, '-k')
    check_option(parser, options.sam_file, '-s')

    if None != options.seed:
        random.seed(options.seed)

    flow_order = options.flow_order.upper()
    key_sequence = options.key_sequence.upper()

    rg_set = 0
    fp = open(options.sam_file, 'r')
    for line in fp:
        line = line.rstrip('\r\n')
        if line[0] == '@':
            if rg_set < 0:
                raise( 'RG was already set' )
            if line[0:3] == '@RG':
                line = line + ('\tFO:%s\tKS:%s' % (flow_order, options.key_sequence)) 
                rg_set = 1
            sys.stdout.write(line + '\n')
            continue
        elif 0 == rg_set:
            sys.stdout.write('@RG\tID:XYZ\tSM:SM\tFO:%s\tKS:%s\n' % (flow_order, options.key_sequence))
        if 0 == rg_set:
            rg_set = -1 # dummy
        elif 0 < rg_set:
            rg_set = -2 # add to previous
        if -1 == rg_set:
            sys.stdout.write(line + '\tRG:Z:XYZ\tFZ:B:S')
        else:
            sys.stdout.write(line + '\tFZ:B:S')
        tokens = line.split('\t')
        seq = tokens[9].upper()
        if tokens[1] == '16':
            seq = reverse_compliment(seq)
        elif tokens[1] != '4' and tokens[1] != '0':
            raise( 'unknown flag' )
        seq = key_sequence + seq
        i = j = 0;
        while i < len(seq):
            while seq[i] != flow_order[j]:
                sys.stdout.write(',' + str(make_flow_signal(0)))
                j = (j + 1) % len(flow_order)
            l = 1
            while i < len(seq)-1 and seq[i] == seq[i+1]:
                l += 1
                i += 1
            sys.stdout.write(',' + str(make_flow_signal(l)))
            j = (j + 1) % len(flow_order)
            i += 1
        sys.stdout.write('\n')
    fp.close()

if __name__ == "__main__":
    __main__()
