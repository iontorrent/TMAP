#!/usr/bin/env python

# TODO
# - map1/map2/map3 specific options

"""
Runs TMAP on Ion Torrent data.
Produces a SAM file containing the mappings.
Works with TMAP version 0.0.8 or higher.

usage: tmap_wrapper.py [options]
    --threads: The number of threads to use
    --ref: The reference genome to use or index
    --fastq: The fastq file to use for the mapping
    --output: The file to save the output (SAM format)
    --params: Parameter setting to use (pre_set or full)
    --fileSource: Whether to use a previously indexed reference sequence or one from history (indexed or history)
    --matchScore: The match score
    --mismatchPenalty: The ismatch penalty
    --gapOpenPenalty: The gap open penalty
    --gapExtensPenalty: The gap extension penalty
    --flowPenalty: The flow score penalty
    --flowOrder: The flow order ([ACGT]{4+} or "sff")
    --bandWidth: The band width
    --globalMap: The soft-clipping type (0 - allow on the right and left, 1 - allow on the left, 2 - allow on the right, 3 - do not allow soft-clipping)
    --duplicateWindow: Remove duplicate alignments from different algorithms within this bp window (-1 to disable)
    --scoringThreshold: The score threshold divided by the match score 
    --queueSize: The queue size for the reads
    --outputFilter: The output filter (0 - unique best hits, 1 - random best hit, 2 - all best htis, 3 - all alignments)
    --rgTag: The flag to specify RG tag(s) in the SAM header
    --rgTagID: The RG ID tag to add to the SAM header
    --rgTagCN: The RG CN tag to add to the SAM header
    --rgTagDS: The RG DS tag to add to the SAM header
    --rgTagDT: The RG DT tag to add to the SAM header
    --rgTagLB: The RG LB tag to add to the SAM header
    --rgTagPI: The RG PI tag to add to the SAM header
    --rgTagPL: The RG PL tag to add to the SAM header
    --rgTagPU: The RG PU tag to add to the SAM header
    --rgTagSM: The RG SM tag to add to the SAM header
    --filterIndependently: Apply the output filter for each algorithm independently

    --map1: Flag to run map1 in the first stage
    --map1SeedLength: The k-mer length to seed CALs (-1 to disable)
    --map1SeedMismatches: The maximum number of mismatches in the seed 
    --map1SecondarySeedLength: The secondary seed length (-1 to disable)
    --map1NumEdits: The maximum number of edits or false-negative probability assuming the maximum error rate
    --map1BaseError: The assumed per-base maximum error rate
    --map1Mismatches: The maximum number of or (read length) fraction of mismatches 
    --map1GapOpens: The maximum number of or (read length) fraction of indel starts 
    --map1GapExtensions: The maximum number of or (read length) fraction of indel extensions 
    --map1MaxCALsDeletion: The maximum number of CALs to extend a deletion 
    --map1EndIndels: Indels are not allowed within this number of bps from the end of the read 
    --map1MaxOptimalCALs: Stop searching when INT optimal CALs have been found 
    --map1MaxNodes: The maximum number of alignment nodes 

    --map2: Flag to run map2 in the first stage
    --map2Coefficient: The coefficient of length-threshold adjustment 
    --map2SeedIntervalSize: The maximum seeding interval size 
    --map2ZBest: Keep the z-best nodes during prefix trie traversal
    --map2ReverseTrigger: The # seeds to trigger reverse alignment 

    --map3: Flag to run map3 in the first stage
    --map3SeedLength: The k-mer length to seed CALs (-1 tunes to the genome size) 
    --map3SeedMaxHits: The maximum number of hits returned by a seed 
    --map3SeedWindow: The window of bases in which to group seeds 
    --map3HPEnumeration: The single homopolymer error difference for enumeration 
    
    --MAP1: Flag to run MAP1 in the second stage
    --MAP1SeedLength: The k-mer length to seed CALs (-1 to disable)
    --MAP1SeedMismatches: The maximum number of mismatches in the seed 
    --MAP1SecondarySeedLength: The secondary seed length (-1 to disable)
    --MAP1NumEdits: The maximum number of edits or false-negative probability assuming the maximum error rate
    --MAP1BaseError: The assumed per-base maximum error rate
    --MAP1Mismatches: The maximum number of or (read length) fraction of mismatches 
    --MAP1GapOpens: The maximum number of or (read length) fraction of indel starts 
    --MAP1GapExtensions: The maximum number of or (read length) fraction of indel extensions 
    --MAP1MaxCALsDeletion: The maximum number of CALs to extend a deletion 
    --MAP1EndIndels: Indels are not allowed within this number of bps from the end of the read 
    --MAP1MaxOptimalCALs: Stop searching when INT optimal CALs have been found 
    --MAP1MaxNodes: The maximum number of alignment nodes 

    --MAP2: Flag to run MAP2 in the second stage
    --MAP2Coefficient: The coefficient of length-threshold adjustment 
    --MAP2SeedIntervalSize: The maximum seeding interval size 
    --MAP2ZBest: Keep the z-best nodes during prefix trie traversal
    --MAP2ReverseTrigger: The # seeds to trigger reverse alignment 

    --MAP3: Flag to run MAP3 in the second stage
    --MAP3SeedLength: The k-mer length to seed CALs (-1 tunes to the genome size) 
    --MAP3SeedMaxHits: The maximum number of hits returned by a seed 
    --MAP3SeedWindow: The window of bases in which to group seeds 
    --MAP3HPEnumeration: The single homopolymer error difference for enumeration 

    --suppressHeader: Suppress header
    --dbkey: Dbkey for reference genome
    --do_not_build_index: Flag to specify that provided file is already indexed and to just use 'as is'
"""

import optparse, os, shutil, subprocess, sys, tempfile

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def map1_parse( map1SeedLength, map1SeedMismatches, \
        map1SecondarySeedLength, map1NumEdits, \
        map1BaseError, map1Mismatches, \
        map1GapOpens, map1GapExtensions, map1MaxCALsDeletion, \
        map1EndIndels, map1MaxOptimalCALs, map1MaxNodes ):
    return '-l %s -s %s -L %s -p %s -P %s -m %s -o %s -e %s -d %s -i %s -b %s -Q %s' % \
            (map1SeedLength, map1SeedMismatches, \
            map1SecondarySeedLength, map1NumEdits, \
            map1BaseError, map1Mismatches, \
            map1GapOpens, map1GapExtensions, map1MaxCALsDeletion, \
            map1EndIndels, map1MaxOptimalCALs, map1MaxNodes )

def map2_parse( map2Coefficient, map2SeedIntervalSize, map2ZBest, map2ReverseTrigger ):
    return '-c %s -S %s -b %s -N %s' % ( map2Coefficient, map2SeedIntervalSize, map2ZBest, map2ReverseTrigger )

def map3_parse( map3SeedLength, map3SeedMaxHits, map3SeedWindow, map3HPEnumeration ):
    return '-l %s -S %s -b %s -H %s' % ( map3SeedLength, map3SeedMaxHits, map3SeedWindow, map3HPEnumeration )

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    # Global options
    parser.add_option( '--threads', dest='threads', help='The number of threads to use' )
    parser.add_option( '--ref', dest='ref', help='The reference genome to use or index' )
    parser.add_option( '--fastq', dest='fastq', help='The fastq file to use for the mapping' )
    parser.add_option( '--output', dest='output', help='The file to save the output (SAM format)' )
    parser.add_option( '--params', dest='params', help='Parameter setting to use (pre_set or full)' )
    parser.add_option( '--fileSource', dest='fileSource', help='Whether to use a previously indexed reference sequence or one from history (indexed or history)' )
    parser.add_option( '--mapall', dest='mapall', help='Flag indicating if mapall options are present')
    parser.add_option( '--matchScore', dest='matchScore', help='The match score' )
    parser.add_option( '--mismatchPenalty', dest='mismatchPenalty', help='Mismatch penalty' )
    parser.add_option( '--gapOpenPenalty', dest='gapOpenPenalty', help='Gap open penalty' )
    parser.add_option( '--gapExtensPenalty', dest='gapExtensPenalty', help='Gap extension penalty' )
    parser.add_option( '--flowPenalty', dest='flowPenalty', help='Flow score penalty' )
    parser.add_option( '--flowOrder', dest='flowOrder', help='Flow order' )
    parser.add_option( '--bandWidth', dest='bandWidth', help='The band width' )
    parser.add_option( '--globalMap', dest='globalMap', help='Map the full read (no soft-clipping)' )
    parser.add_option( '--duplicateWindow', dest='duplicateWindow', help='Remove duplicate alignments from different algorithms within this bp window (-1 to disable)' )
    parser.add_option( '--scoringThreshold', dest='scoringThreshold', help='The score threshold divided by the match score ' )
    parser.add_option( '--queueSize', dest='queueSize', help='The queue size for the reads' )
    parser.add_option( '--outputFilter', dest='outputFilter', help='The output filter (0 - unique best hits, 1 - random best hit, 2 - all best htis, 3 - all alignments)' )
    parser.add_option( '--rgTag', dest='rgTag', help='The flag to specify RG tag(s) in the SAM header' )
    parser.add_option( '--rgTagID', dest='rgTagID', default='', help='The RG ID to add to the SAM header' )
    parser.add_option( '--rgTagCN', dest='rgTagCN', default='', help='The RG CN to add to the SAM header' )
    parser.add_option( '--rgTagDS', dest='rgTagDS', default='', help='The RG DS to add to the SAM header' )
    parser.add_option( '--rgTagDT', dest='rgTagDT', default='', help='The RG DT to add to the SAM header' )
    parser.add_option( '--rgTagLB', dest='rgTagLB', default='', help='The RG LB to add to the SAM header' )
    parser.add_option( '--rgTagPI', dest='rgTagPI', default='', help='The RG PI to add to the SAM header' )
    parser.add_option( '--rgTagPL', dest='rgTagPL', default='', help='The RG PL to add to the SAM header' )
    parser.add_option( '--rgTagPU', dest='rgTagPU', default='', help='The RG PU to add to the SAM header' )
    parser.add_option( '--rgTagSM', dest='rgTagSM', default='', help='The RG SM to add to the SAM header' )
    parser.add_option( '--filterIndependently', dest='filterIndependently', help='Apply the output filter for each algorithm independently' )
    parser.add_option( '--suppressHeader', dest='suppressHeader', help='Suppress header' )
    parser.add_option( '--dbkey', dest='dbkey', help='Dbkey for reference genome' )
    parser.add_option( '--do_not_build_index', dest='do_not_build_index', action='store_true', help="Don't build index" )
    # map 1 - stage 1
    parser.add_option( '--map1', dest='map1', help='True if map1 should be run in the first stage' )
    parser.add_option( '--map1SeedLength', dest='map1SeedLength', help='The k-mer length to seed CALs (-1 to disable)' )
    parser.add_option( '--map1SeedMismatches', dest='map1SeedMismatches', help='The maximum number of mismatches in the seed ' )
    parser.add_option( '--map1SecondarySeedLength', dest='map1SecondarySeedLength', help='The secondary seed length (-1 to disable)' )
    parser.add_option( '--map1NumEdits', dest='map1NumEdits', help='The maximum number of edits or false-negative probability assuming the maximum error rate' )
    parser.add_option( '--map1BaseError', dest='map1BaseError', help='The assumed per-base maximum error rate' )
    parser.add_option( '--map1Mismatches', dest='map1Mismatches', help='The maximum number of or (read length) fraction of mismatches ' )
    parser.add_option( '--map1GapOpens', dest='map1GapOpens', help='The maximum number of or (read length) fraction of indel starts ' )
    parser.add_option( '--map1GapExtensions', dest='map1GapExtensions', help='The maximum number of or (read length) fraction of indel extensions ' )
    parser.add_option( '--map1MaxCALsDeletion', dest='map1MaxCALsDeletion', help='The maximum number of CALs to extend a deletion ' )
    parser.add_option( '--map1EndIndels', dest='map1EndIndels', help='Indels are not allowed within this number of bps from the end of the read ' )
    parser.add_option( '--map1MaxOptimalCALs', dest='map1MaxOptimalCALs', help='Stop searching when INT optimal CALs have been found ' )
    parser.add_option( '--map1MaxNodes', dest='map1MaxNodes', help='The maximum number of alignment nodes ' )
    # map 2 - stage 1
    parser.add_option( '--map2', dest='map2', help='True if map2 should be run in the first stage' )
    parser.add_option( '--map2Coefficient', dest='map2Coefficient', help='The coefficient of length-threshold adjustment ' )
    parser.add_option( '--map2SeedIntervalSize', dest='map2SeedIntervalSize', help='The maximum seeding interval size ' )
    parser.add_option( '--map2ZBest', dest='map2ZBest', help='Keep the z-best nodes during prefix trie traversal' )
    parser.add_option( '--map2ReverseTrigger', dest='map2ReverseTrigger', help='The # seeds to trigger reverse alignment ' )
    # map 3 - stage 1
    parser.add_option( '--map3', dest='map3', help='True if map3 should be run in the first stage' )
    parser.add_option( '--map3SeedLength', dest='map3SeedLength', help='The k-mer length to seed CALs (-1 tunes to the genome size) ' )
    parser.add_option( '--map3SeedMaxHits', dest='map3SeedMaxHits', help='The maximum number of hits returned by a seed ' )
    parser.add_option( '--map3SeedWindow', dest='map3SeedWindow', help='The window of bases in which to group seeds ' )
    parser.add_option( '--map3HPEnumeration', dest='map3HPEnumeration', help='The single homopolymer error difference for enumeration ' )
    # map 1 - stage 2
    parser.add_option( '--MAP1', dest='MAP1', help='True if map1 should be run in the second stage' )
    parser.add_option( '--MAP1SeedLength', dest='MAP1SeedLength', help='The k-mer length to seed CALs (-1 to disable)' )
    parser.add_option( '--MAP1SeedMismatches', dest='MAP1SeedMismatches', help='The maximum number of mismatches in the seed ' )
    parser.add_option( '--MAP1SecondarySeedLength', dest='MAP1SecondarySeedLength', help='The secondary seed length (-1 to disable)' )
    parser.add_option( '--MAP1NumEdits', dest='MAP1NumEdits', help='The maximum number of edits or false-negative probability assuming the maximum error rate' )
    parser.add_option( '--MAP1BaseError', dest='MAP1BaseError', help='The assumed per-base maximum error rate' )
    parser.add_option( '--MAP1Mismatches', dest='MAP1Mismatches', help='The maximum number of or (read length) fraction of mismatches ' )
    parser.add_option( '--MAP1GapOpens', dest='MAP1GapOpens', help='The maximum number of or (read length) fraction of indel starts ' )
    parser.add_option( '--MAP1GapExtensions', dest='MAP1GapExtensions', help='The maximum number of or (read length) fraction of indel extensions ' )
    parser.add_option( '--MAP1MaxCALsDeletion', dest='MAP1MaxCALsDeletion', help='The maximum number of CALs to extend a deletion ' )
    parser.add_option( '--MAP1EndIndels', dest='MAP1EndIndels', help='Indels are not allowed within this number of bps from the end of the read ' )
    parser.add_option( '--MAP1MaxOptimalCALs', dest='MAP1MaxOptimalCALs', help='Stop searching when INT optimal CALs have been found ' )
    parser.add_option( '--MAP1MaxNodes', dest='MAP1MaxNodes', help='The maximum number of alignment nodes ' )
    # map 2 - stage 2
    parser.add_option( '--MAP2', dest='MAP2', help='True if map2 should be run in the second stage' )
    parser.add_option( '--MAP2Coefficient', dest='MAP2Coefficient', help='The coefficient of length-threshold adjustment ' )
    parser.add_option( '--MAP2SeedIntervalSize', dest='MAP2SeedIntervalSize', help='The maximum seeding interval size ' )
    parser.add_option( '--MAP2ZBest', dest='MAP2ZBest', help='Keep the z-best nodes during prefix trie traversal' )
    parser.add_option( '--MAP2ReverseTrigger', dest='MAP2ReverseTrigger', help='The # seeds to trigger reverse alignment ' )
    # map 3 - stage 2
    parser.add_option( '--MAP3', dest='MAP3', help='True if map3 should be run in the second stage' )
    parser.add_option( '--MAP3SeedLength', dest='MAP3SeedLength', help='The k-mer length to seed CALs (-1 tunes to the genome size) ' )
    parser.add_option( '--MAP3SeedMaxHits', dest='MAP3SeedMaxHits', help='The maximum number of hits returned by a seed ' )
    parser.add_option( '--MAP3SeedWindow', dest='MAP3SeedWindow', help='The window of bases in which to group seeds ' )
    parser.add_option( '--MAP3HPEnumeration', dest='MAP3HPEnumeration', help='The single homopolymer error difference for enumeration ' )
    # parse the options
    (options, args) = parser.parse_args()

    # output version # of tool
    try:
        tmp = tempfile.NamedTemporaryFile().name
        tmp_stdout = open( tmp, 'wb' )
        proc = subprocess.Popen( args='tmap --version 2>&1', shell=True, stdout=tmp_stdout )
        tmp_stdout.close()
        returncode = proc.wait()
        stdout = None
        for line in open( tmp_stdout.name, 'rb' ):
            if line.lower().find( 'version' ) >= 0:
                stdout = line.strip()
                break
        if stdout:
            sys.stdout.write( 'TMAP %s\n' % stdout )
        else:
            raise Exception
    except:
        sys.stdout.write( 'Could not determine TMAP version\n' )

    fastq = options.fastq

    # make temp directory for placement of indices
    tmp_index_dir = tempfile.mkdtemp()
    tmp_dir = tempfile.mkdtemp()

    # index if necessary
    if options.fileSource == 'history' and not options.do_not_build_index:
        ref_file = tempfile.NamedTemporaryFile( dir=tmp_index_dir )
        ref_file_name = ref_file.name
        ref_file.close()
        os.symlink( options.ref, ref_file_name )
        cmd1 = 'tmap index -f %s -v ' % ( ref_file_name )
        try:
            tmp = tempfile.NamedTemporaryFile( dir=tmp_index_dir ).name
            tmp_stderr = open( tmp, 'wb' )
            proc = subprocess.Popen( args=cmd1, shell=True, cwd=tmp_index_dir, stderr=tmp_stderr.fileno() )
            returncode = proc.wait()
            tmp_stderr.close()
            # get stderr, allowing for case where it's very large
            tmp_stderr = open( tmp, 'rb' )
            stderr = ''
            buffsize = 1048576
            try:
                while True:
                    stderr += tmp_stderr.read( buffsize )
                    if not stderr or len( stderr ) % buffsize != 0:
                        break
            except OverflowError:
                pass
            tmp_stderr.close()
            if returncode != 0:
                raise Exception, stderr
        except Exception, e:
            # clean up temp dirs
            if os.path.exists( tmp_index_dir ):
                shutil.rmtree( tmp_index_dir )
            if os.path.exists( tmp_dir ):
                shutil.rmtree( tmp_dir )
            stop_err( 'Error indexing reference sequence. ' + str( e ) )
    else:
        ref_file_name = options.ref

    # set up mapping and generate mapping command options
    if options.params == 'pre_set':
        mapall_options = '-n %s' % ( options.threads )
        map1_options = 'map1'
        map2_options = 'map2'
        map3_options = 'map3'
        MAP1_options = ''
        MAP2_options = ''
        MAP3_options = ''
    else:
        # mapall options
        if options.rgTag == 'true':
            rgTag = ''
            if options.rgTagID != '':
                rgTag += '-R "ID:' + options.rgTagID + '" '
            if options.rgTagCN != '':
                rgTag += '-R "CN:' + options.rgTagCN + '" '
            if options.rgTagDS != '':
                rgTag += '-R "DS:' + options.rgTagDS + '" '
            if options.rgTagDT != '':
                rgTag += '-R "DT:' + options.rgTagDT + '" '
            if options.rgTagLB != '':
                rgTag += '-R "LB:' + options.rgTagLB + '" '
            if options.rgTagPI != '':
                rgTag += '-R "PI:' + options.rgTagPI + '" '
            if options.rgTagPL != '':
                rgTag += '-R "PL:' + options.rgTagPL + '" '
            if options.rgTagPU != '':
                rgTag += '-R "PU:' + options.rgTagPU + '" '
            if options.rgTagSM != '':
                rgTag += '-R "SM:' + options.rgTagSM + '" '
            rgTag.rstrip(' ')
        else:
            rgTag = ''
        if None != options.flowOrder and '' != options.flowOrder:
            flowOrder = '-x ' + options.flowOrder
        else:
            flowOrder = ''
        if options.filterIndependently == 'true':
            filterIndependently = '-I'
        else:
            filterIndependently = ''
        if options.mapall == 'true':
            mapall_options = '-A %s -M %s -O %s -E %s -X %s %s %s -W %s -T %s -q %s -n %s -a %s %s %s' % \
                    ( options.matchScore, options.mismatchPenalty, options.gapOpenPenalty, options.gapExtensPenalty,
                            options.flowPenalty, flowOrder,
                            options.globalMap, options.duplicateWindow, options.scoringThreshold, options.queueSize,
                            options.threads, options.outputFilter, rgTag, filterIndependently )
        else:
            mapall_options = ''

        # map1 - stage one
        if options.map1 == 'true':
            map1_options = 'map1 %s' % \
                    map1_parse( options.map1SeedLength, options.map1SeedMismatches, \
                    options.map1SecondarySeedLength, options.map1NumEdits, \
                    options.map1BaseError, options.map1Mismatches, \
                    options.map1GapOpens, options.map1GapExtensions, options.map1MaxCALsDeletion, \
                    options.map1EndIndels, options.map1MaxOptimalCALs, options.map1MaxNodes )
        else:
            map1_options = ''
        # map2 - stage one
        if options.map2 == 'true':
            map2_options = 'map2 %s' % \
                    map2_parse( options.map2Coefficient, options.map2SeedIntervalSize, \
                    options.map2ZBest, options.map2ReverseTrigger )
        else:
            map2_options = ''
        # map3 - stage one
        if options.map3 == 'true':
            map3_options = 'map3 %s' % \
                    map3_parse( options.map3SeedLength, options.map3SeedMaxHits, \
                    options.map3SeedWindow, options.map3HPEnumeration )
        else:
            map3_options = ''
        # map1 - stage two
        if options.MAP1== 'true':
            MAP1_options = 'MAP1 %s' % \
                    map1_parse( options.MAP1SeedLength, options.MAP1SeedMismatches, \
                    options.MAP1SecondarySeedLength, options.MAP1NumEdits, \
                    options.MAP1BaseError, options.MAP1Mismatches, \
                    options.MAP1GapOpens, options.MAP1GapExtensions, options.MAP1MaxCALsDeletion, \
                    options.MAP1EndIndels, options.MAP1MaxOptimalCALs, options.MAP1MaxNodes )
        else:
            MAP1_options = ''
        # map2 - stage two
        if options.MAP2 == 'true':
            MAP2_options = 'MAP2 %s' % \
                    map2_parse( options.MAP2Coefficient, options.MAP2SeedIntervalSize, \
                    options.MAP2ZBest, options.MAP2ReverseTrigger )
        else:
            MAP2_options = ''
        # map3 - stage two
        if options.MAP3 == 'true':
            MAP3_options = 'MAP3 %s' % \
                    MAP3_parse( options.MAP3SeedLength, options.MAP3SeedMaxHits, \
                    options.MAP3SeedWindow, options.MAP3HPEnumeration )
        else:
            MAP3_options = ''

    #mapping_cmds 
    # prepare actual mapping and generate mapping commands
    cmd2 = 'tmap mapall -f %s -r %s -F fastq %s -v %s %s %s %s %s %s > %s' % \
            ( ref_file_name, fastq, mapall_options, 
                    map1_options, map2_options, map3_options,
                    MAP1_options, MAP2_options, MAP3_options,
                    options.output )
    # perform alignments
    buffsize = 1048576
    try:
        # need to nest try-except in try-finally to handle 2.4
        try:
            # align
            try:
                tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
                tmp_stderr = open( tmp, 'wb' )
                proc = subprocess.Popen( args=cmd2, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
                returncode = proc.wait()
                tmp_stderr.close()
                # get stderr, allowing for case where it's very large
                tmp_stderr = open( tmp, 'rb' )
                stderr = ''
                try:
                    while True:
                        stderr += tmp_stderr.read( buffsize )
                        if not stderr or len( stderr ) % buffsize != 0:
                            break
                except OverflowError:
                    pass
                tmp_stderr.close()
                if returncode != 0:
                    raise Exception, stderr
            except Exception, e:
                raise Exception, 'Error mapping sequence. ' + str( e )
            # remove header if necessary
            if options.suppressHeader == 'true':
                tmp_out = tempfile.NamedTemporaryFile( dir=tmp_dir)
                tmp_out_name = tmp_out.name
                tmp_out.close()
                try:
                    shutil.move( options.output, tmp_out_name )
                except Exception, e:
                    raise Exception, 'Error moving output file before removing headers. ' + str( e )
                fout = file( options.output, 'w' )
                for line in file( tmp_out.name, 'r' ):
                    if not ( line.startswith( '@HD' ) or line.startswith( '@SQ' ) or line.startswith( '@RG' ) or line.startswith( '@PG' ) or line.startswith( '@CO' ) ):
                        fout.write( line )
                fout.close()
            # check that there are results in the output file
            if os.path.getsize( options.output ) > 0:
                sys.stdout.write( 'TMAP completed' )
            else:
                raise Exception, 'The output file is empty. You may simply have no matches, or there may be an error with your input file or settings.'
        except Exception, e:
            stop_err( 'The alignment failed.\n' + str( e ) )
    finally:
        # clean up temp dir
        if os.path.exists( tmp_index_dir ):
            shutil.rmtree( tmp_index_dir )
        if os.path.exists( tmp_dir ):
            shutil.rmtree( tmp_dir )

if __name__=="__main__": __main__()
