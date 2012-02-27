#!/usr/bin/env python

# TODO
# - map1/map2/map3 specific options

"""
Runs TMAP on Ion Torrent data.
Produces a SAM file containing the mappings.
Works with TMAP version 0.3.3 or higher.

usage: tmap_wrapper.py [options]
    --threads: The number of threads to use
    --ref: The reference genome to use or index
    --input: The input FASTQ/SFF file to use for the mapping
    --inputtype: The input type (FASTQ/SFF)
    --output: The file to save the output (SAM format)
    --params: Parameter setting to use (pre_set or full)
    --fileSource: Whether to use a previously indexed reference sequence or one from history (indexed or history)
    --algorithm: The algorithm (ex. mapall, map1, map2, map3, map4, ...)
    --globalOptions: The global options 
    --flowspaceOptions: The flowspace options
    --pairingOptions: The pairing options
    --algorithmOptions: The algorithm options
    --suppressHeader: Suppress header
    --dbkey: Dbkey for reference genome
    --do_not_build_index: Flag to specify that provided file is already indexed and to just use 'as is'
"""

import optparse, os, shutil, subprocess, sys, tempfile

def stop_err( msg ):
    sys.stderr.write( '%s\n' % msg )
    sys.exit()

def __main__():
    #Parse Command Line
    parser = optparse.OptionParser()
    # Global options
    parser.add_option( '--threads', dest='threads', help='The number of threads to use' )
    parser.add_option( '--ref', dest='ref', help='The reference genome to use or index' )
    parser.add_option( '--input', dest='input', help='The input file to use for the mapping' )
    parser.add_option( '--inputtype', dest='inputtype', help='The input file type' )
    parser.add_option( '--output', dest='output', help='The file to save the output (SAM format)' )
    parser.add_option( '--params', dest='params', help='Parameter setting to use (pre_set or full)' )
    parser.add_option( '--fileSource', dest='fileSource', help='Whether to use a previously indexed reference sequence or one from history (indexed or history)' )
    parser.add_option( '--algorithm', dest='algorithm', help='The algorithm (ex. mapall, map1, map2, map3, map4, ...)')
    parser.add_option( '--globalOptions', dest='globalOptions', help='The global options ' )
    parser.add_option( '--flowspaceOptions', dest='flowspaceOptions', help='The flowspace options' )
    parser.add_option( '--pairingOptions', dest='pairingOptions', help='The pairing options' )
    parser.add_option( '--algorithmOptions', dest='algorithmOptions', help='The algorithm options' )
    parser.add_option( '--suppressHeader', dest='suppressHeader', help='Suppress header' )
    parser.add_option( '--dbkey', dest='dbkey', help='Dbkey for reference genome' )
    parser.add_option( '--do_not_build_index', dest='do_not_build_index', action='store_true', help="Don't build index" )
    
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
        options.algorithm = 'mapall'
        options.flowspaceOptions = ''
        options.pairingOptions = ''
        options.algorithmOptions = 'stage1 map1 map2 map3'

    #mapping_cmds 
    # prepare actual mapping and generate mapping commands
    cmd = 'tmap %s -f %s -r %s -i %s -n %s %s %s %s' % \
            ( options.algorithm, options.ref, options.input, options.inputtype, \
            options.threads, options.flowspaceOptions, \
            options.pairingOptions, options.algorithmOptions ) 
    # perform alignments
    buffsize = 1048576
    try:
        # need to nest try-except in try-finally to handle 2.4
        try:
            # align
            try:
                tmp = tempfile.NamedTemporaryFile( dir=tmp_dir ).name
                tmp_stderr = open( tmp, 'wb' )
                proc = subprocess.Popen( args=cmd, shell=True, cwd=tmp_dir, stderr=tmp_stderr.fileno() )
                returncode = proc.wait()
                tmp_stderr.close()
                # get stderr, allowing for case where it's very large
                tmp_stderr = open( tmp, 'rb' )
                stderr = ''
                try:
                    stderr += cmd + '\n'
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
