#!/usr/bin/env python
import sys
import os
import subprocess

import reformat_contacts
from localconfig import *


def check_output(command, wfile, ok_to_fail=False):
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if ok_to_fail or not p.returncode:
        return stdout
    msg = "ERROR: Command %s exited with nonzero returncode %s. Stdout available in %s, stderr in %s"
    err = wfile + '.fail.stderr'
    out = wfile + '.fail.stdout'
    print >> open(err, 'w'), stderr
    print >> open(out, 'w'), stdout
    print >> sys.stderr, msg % (command, p.returncode, out, err)
    sys.exit(1)


def main(seqfile, cfile, factor=1.0, nohoms_flag=False):
    saved_path = os.getcwd()
    rundir = os.path.dirname(os.path.abspath(seqfile)) + '/'

    if not os.path.exists(rundir + 'rosetta'):
        check_output(['mkdir', rundir + 'rosetta'], seqfile)

    ### reformat the contact maps into constraint files used by rosetta
    ### creates: "cfile-<factor>.constraints"
    ros_cfile = cfile + '-' + str(factor) + '.constraints'

    reformat_contacts.reformat(seqfile, cfile, factor, outfile_name=ros_cfile)

    ### run make_fragments.pl to create fragment database
    check_output(['cp', seqfile, rundir + 'rosetta'], seqfile)
    os.chdir(rundir + 'rosetta')
    if nohoms_flag:
        check_output([rosetta_make_fragments, os.path.basename(seqfile), '-nosam', '-noporter', '-nohoms'], seqfile)
    else:
        check_output([rosetta_make_fragments, os.path.basename(seqfile), '-nosam', '-noporter'], seqfile)

    os.chdir(saved_path)


if __name__ == "__main__":
    
    nohoms_flag = False
    factor = 1.0

    if '-f' in sys.argv:
        idx = sys.argv.index('-f')
        try:
            factor = int(sys.argv[idx+1])
        except:
            print 'Factor of sequence length (determining number of constraints to be used during folding) -f must be float, %r is not. Default is %s.' % (sys.argv[idx+1], factor)
            sys.exit(1)
        del sys.argv[idx]
        del sys.argv[idx]

    if '--nohoms' in sys.argv:
        idx = sys.argv.index('--nohoms')
        nohoms_flag = True
        del sys.argv[idx]

    seqfile = sys.argv[1]
    cfile = sys.argv[2]

    main(seqfile, cfile, factor=factor, nohoms_flag=nohoms_flag)
