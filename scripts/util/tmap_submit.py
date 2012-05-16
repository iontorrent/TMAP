#!/usr/bin/env python
# Copyright (C) 2010 Ion Torrent Systems, Inc. All Rights Reserved

import sys
import os.path
import re
from optparse import OptionParser
from subprocess import call

def check_option(parser, value, name):
    if None == value:
        print 'Option ' + name + ' required.\n'
        parser.print_help()
        sys.exit(1)

def check_file(file_name):
    if not os.path.isfile(file_name):
        print 'File not found: ' + str(file_name)
        sys.exit(1)

def check_dir(dir_name):
    if not os.path.isdir(dir_name):
        print 'Directory not found: ' + str(dir_name)
        print 'making directory: ' + str(dir_name)
        try:
            os.mkdir(dir_name)
        except:
            print "Directory creation failed.  You apparently don't have sufficient write permissions to: " + str(dir_name)
            sys.exit(1)

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('-f', '--fn-fasta', help='the FASTA reference file name', dest='fn_fasta')
    parser.add_option('-r', '--fn-reads', help='the reads file name', dest='fn_reads')
    parser.add_option('-s', '--sam-directory', help='where the SAM file should be made', dest='sam_dir')
    parser.add_option('-n', '--num-threads', help='the number of threads', dest='num_threads')
    parser.add_option('--tmap-path', help='the path to the tmap binary', dest='tmap_path')
    parser.add_option('--dwgsim-eval-path', help='the path to the dwgsim_eval binary (if the reads were simulated with dwgsim_eval)', dest='dwgsim_eval_path', default='')
    parser.add_option('--alignStats-path', help='the path to the alignStast binary', dest='alignStats_path', default='')
    parser.add_option('--mapping-algorithm', help='the mapping command (ex. mapall, map1, map2 map3)', dest='mapping_algorithm')
    parser.add_option('--mapall-algorithms', help='the mapping algorithms and options for mapall (only used when --mapping-algorithm=mapall)', dest='mapall_algorithms')
    parser.add_option('--sam-sff-tags', help='include SFF specific SAM tags', action="store_true", default=False, dest='sam_sff_tags')
    parser.add_option('--pbs-queue', help='the PBS queue to which to submit')
    parser.add_option('--redirect-stdout', help='redirects stdout of tmap to the specified file', dest="tmap_stdout", default='')
    parser.add_option('--profile-tmap', help='profiles tmap and requires google-perftools be in your LD_LIBRARY_PATH', default='', dest="profile_tmap")
    parser.add_option('--submit-script-dir', help="name for submission script, override default", dest="submit_script_dir")
    parser.add_option('--additional-pbs-settings', help="comman seperated list of additioanl PBS/environemnt settings.  formating is up to you, user.", dest="pbs_settings", default='')
    (options, args) = parser.parse_args()
    if len(args) != 0:
        parser.print_help()
        sys.exit(1)

    # Check options
    check_option(parser, options.fn_fasta, '-f')
    check_option(parser, options.fn_reads, '-r')
    check_option(parser, options.sam_dir, '-s')
    check_option(parser, options.num_threads, '-n')
    check_option(parser, options.tmap_path, '--tmap-path')
    check_option(parser, options.dwgsim_eval_path, '--dwgsim_eval_path')
    if '' == options.dwgsim_eval_path:
        print "Warning: dwgsim_eval will not be run."
    check_option(parser, options.alignStats_path, '--alignStats-path')
    if '' == options.alignStats_path:
        print "Warning: alignStats will not be run."
    check_option(parser, options.mapping_algorithm, '--mapping-algorithm')
    if 'mapall' == options.mapping_algorithm.split()[0]:
        mapall_use = True
        check_option(parser, options.mapall_algorithms, '--mapall-algorithms')
    else:
        mapall_use = False
        options.mapall_algorithms = ''
    check_option(parser, options.pbs_queue, '--pbs-queue')
    
    if options.tmap_stdout != '':
        options.tmap_stdout = (" > %s" ) % (options.tmap_stdout)
    # Genome info file
    p = re.compile('.fasta$')
    fn_fasta_info = p.sub('.info.txt', options.fn_fasta)

    # Check paths
    check_file(options.fn_fasta)
    check_file(fn_fasta_info)
    check_file(options.fn_reads)
    check_dir(options.sam_dir)
    check_file(options.tmap_path)
    if '' != options.dwgsim_eval_path:
        check_file(options.dwgsim_eval_path)
    if '' != options.alignStats_path:
        check_file(options.alignStats_path)
    
    # Check mapping algorithms
    if not options.mapping_algorithm.split()[0] in ('mapall', 'map1', 'map2', 'map3'):
        #just in case it's something like
        #map3 --seed-length 55                                                       
        print 'Unrecognized mapping algorithm: ' + options.mapping_algorithm
    if '' != options.profile_tmap:
        options.profile_tmap = "CPUPROFILE=" + options.profile_tmap
    # Create the tmap command 
    if mapall_use:
        temp_mapall_algorithms = options.mapall_algorithms
        temp = 'mapall_' + ('_'.join(temp_mapall_algorithms.split(' ')))
    else:
        temp_mapping_algorithm = options.mapping_algorithm
        temp = temp_mapping_algorithm.split()[0] + ('_'.join(temp_mapping_algorithm.split(' ')[1:]))
    fn_sam = options.sam_dir + "/" + os.path.basename(options.fn_reads) + "." + temp + ".sam"
    fn_script = options.submit_script_dir + "/" + os.path.basename(options.fn_reads) + "." + temp + ".sh"
    sam_sff_tags = ''
    if options.sam_sff_tags:
        # TODO: check that it is an SFF
        sam_sff_tags = '-Y'
    tmap_cmd = "%s time -p %s %s -f %s -r %s -s %s -n %s -v %s %s %s" % (
            options.profile_tmap,
            options.tmap_path,
            options.mapping_algorithm,
            options.fn_fasta,
            options.fn_reads,
            fn_sam,
            options.num_threads,
            sam_sff_tags,
            options.mapall_algorithms,
            options.tmap_stdout
            )
    # Create the dwgsim_eval command
    dwgsim_eval_cmd = ''
    if '' != options.dwgsim_eval_path:
        fn_eval_txt = fn_sam+ ".eval.txt"
        dwgsim_eval_cmd = "%s -g 10 -z -S %s > %s" % (
                options.dwgsim_eval_path,
                fn_sam,
                fn_eval_txt)

    # Create the alignStats command
    alignStats_cmd = ''
    if '' != options.alignStats_path:
        alignStats_cmd = "%s -a /dev/null -x -i %s -o %s -g %s -q 7,10,17,20,47 -n %s -B 100000" % (
                options.alignStats_path,
                fn_sam,
                fn_sam,
                fn_fasta_info,
                options.num_threads)

    # Create the shell script
    fp_script = open(fn_script, 'w')
    fp_script.write('#!/usr/bin/env bash\n')
    print options.pbs_settings
    if '' != options.pbs_settings:
        #writes all additional pbs settings right here.
        #formatting is left up to the user
        fp_script.write('\n'.join(options.pbs_settings.split(',')) + '\n')
    fp_script.write(
            """
#*! @function
#  @param  $*  the command to be executed
run ()
{
    echo "running: $*";
    eval $*;
    EXIT_CODE="$?";
    if test ${EXIT_CODE} != 0; then
        echo "Error: running '$*'";
        exit 1;
    fi
}
"""
)
    fp_script.write('# TMAP\nrun "%s"\n\n' % tmap_cmd)
    if '' != dwgsim_eval_cmd:
        fp_script.write('# DWGSIM_EVAL\nrun "%s"\n\n' % dwgsim_eval_cmd)
    if '' != alignStats_cmd:
        fp_script.write('# ALIGNSTATS\nrun "%s"\n\n' % alignStats_cmd)
    fp_script.close()

    # Submit the job
    fn_script_stderr = fn_script + '.stderr'
    fn_script_stdout = fn_script + '.stdout'
    returncode = call(['qsub', 
        '-q', options.pbs_queue, 
        '-d', os.getcwd(), 
        '-o', fn_script_stdout, 
        '-e', fn_script_stderr, 
        '-l', 'nodes=1:ppn=' + options.num_threads,
        fn_script
        ])
    if 0 != returncode:
        print "Error submitting the job"
        sys.exit(1)
