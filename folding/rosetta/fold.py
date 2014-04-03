#!/usr/bin/env python
import sys
import os
import math
from subprocess import call
from subprocess import Popen

from localconfig import *


def main(seqfile, ros_cfile, n_cores=1, n_decoys=2000, rundir_postfix='', fragfile_3mers='', fragfile_9mers='', ss2_file=''):

    saved_path = os.getcwd()
    rundir = os.path.dirname(os.path.abspath(seqfile)) + '/'

    ### The number of decoy structures produced in a single run depends on the 
    ### of cores. In total it should be at least 2000 to get native-like structures.
    ### (see PconsFold paper)
    decoys_per_core = int(math.ceil(float(n_decoys) / n_cores))

    os.chdir(rundir + rundir_postfix)

    ### n_cores folders containing run 1..n_cores
    plist = []
    for core in range(1, n_cores + 1):
        optionstr = ''
        optionstr += '-in:file:fasta ' + os.path.abspath(seqfile) + '\n'
        if fragfile_3mers:
            optionstr += '-in:file:frag3 ' + os.path.abspath(fragfile_3mers) + '\n'
        else:
            optionstr += '-in:file:frag3 ' + os.path.abspath('t001_.200.3mers') + '\n'
        if fragfile_9mers:
            optionstr += '-in:file:frag9 ' + os.path.abspath(fragfile_9mers) + '\n'
        else:
            optionstr += '-in:file:frag9 ' + os.path.abspath('t001_.200.9mers') + '\n'
        optionstr += '-abinitio:use_filters true\n'
        optionstr += '-constraints:cst_file ' + os.path.abspath(ros_cfile) + '\n'
        optionstr += '-database ' + rosetta_db_dir + '\n'
        if ss2_file:
            optionstr += '-psipred_ss2 ' + os.path.abspath(ss2_file) + '\n'
        else:
            optionstr += '-psipred_ss2 ' + os.path.abspath('t001_.psipred_ss2') + '\n'
        optionstr += '-out:nstruct ' + str(decoys_per_core) + '\n'
        optionstr += '-seed_offset ' + str(core) + '\n'

        call(['mkdir', 'run_' + str(core)])
        os.chdir('run_' + str(core))
        optionfile = open('abinitio.options', 'w')
        optionfile.write(optionstr)
        optionfile.close()

        plist.append(Popen([rosetta_abinitiorelax, '@abinitio.options']))
        os.chdir('..')

    for p in plist:
        p.wait()

    os.chdir(saved_path)



if __name__ == "__main__":

    ### default value
    n_decoys = 2000
    rundir_postfix = ''
    fragfile_3mers = ''
    fragfile_9mers = ''
    ss2_file = ''

    ### parse parameters
    if '-c' in sys.argv:
        idx = sys.argv.index('-c')    
        try:
            n_cores = int(sys.argv[idx+1])
        except:
            print 'Number of cores -c must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_cores)
            sys.exit(1)
        del sys.argv[idx]
        del sys.argv[idx]

    if '-n' in sys.argv:
        idx = sys.argv.index('-n')    
        try:
            n_decoys = int(sys.argv[idx+1])
        except:
            print 'Number of decoys -n must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_decoys)
            sys.exit(1)
        del sys.argv[idx]
        del sys.argv[idx]

    if '--rundir_postfix' in sys.argv:
        idx = sys.argv.index('--rundir_postfix')    
        rundir_postfix = sys.argv[idx+1]
        del sys.argv[idx]
        del sys.argv[idx]

    if '--3mers' in sys.argv:
        idx = sys.argv.index('--3mers')    
        fragfile_3mers = sys.argv[idx+1]
        del sys.argv[idx]
        del sys.argv[idx]

    if '--9mers' in sys.argv:
        idx = sys.argv.index('--9mers')    
        fragfile_9mers = sys.argv[idx+1]
        del sys.argv[idx]
        del sys.argv[idx]

    if '--psipred_ss2' in sys.argv:
        idx = sys.argv.index('--psipred_ss2')    
        ss2_file = sys.argv[idx+1]
        del sys.argv[idx]
        del sys.argv[idx]

    seqfile = os.path.abspath(sys.argv[1])
    ros_cfile = os.path.abspath(sys.argv[2])

    main(seqfile, ros_cfile, n_cores=n_cores, n_decoys=n_decoys, rundir_postfix=rundir_postfix, fragfile_3mers=fragfile_3mers, fragfile_9mers=fragfile_9mers, ss2_file=ss2_file)
