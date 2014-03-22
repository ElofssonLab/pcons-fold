import sys
import os
from subprocess import call

from reformat_contacts import reformat
from localconfig import *

if __name__ == '__main__':

    os.chdir(rootdir)

    seqfile = sys.argv[1]
    cfile = sys.argv[2]
    factor = sys.argv[3]

    call(['mkdir', rundir])

    ### reformat the contact maps into constraint files used by rosetta
    ### creates: "cfile-<factor>.constraints"
    reformat(seqfile, cfile, factor)
    ros_cfile = cfile[:-4] + '-' + str(factor) + '.constraints'
    call(['mv', ros_cfile, rundir])

    ### run make_fragments.pl to create fragment database
    ### NOTICE: paths in make_fragments.pl need to be adjusted
    call(['cp', seqfile, rundir])
    call(['cp', seqfile + '.ss2', rundir])
    os.chdir(rundir)
    call([rosettadir + '/rosetta_fragments/make_fragments.pl', seqfile, '-nosam', '-noporter', '-psipredfile', seqfile + '.ss2', '-id', '1xxxx'])


