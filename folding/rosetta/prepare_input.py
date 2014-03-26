import sys
import os
from subprocess import call

import reformat_contacts
from localconfig import *

if __name__ == '__main__':

    seqfile = sys.argv[1]
    cfile = sys.argv[2]
    factor = sys.argv[3]

    if not os.path.exists('rosetta'):
        call(['mkdir', 'rosetta'])

    ### reformat the contact maps into constraint files used by rosetta
    ### creates: "cfile-<factor>.constraints"
    ros_cfile = cfile + '-' + str(factor) + '.constraints'
    
    reformat_contacts.reformat(seqfile, cfile, factor, outfile_name=ros_cfile)
    call(['mv', ros_cfile, 'rosetta'])

    ### run make_fragments.pl to create fragment database
    ### NOTICE: paths in make_fragments.pl need to be adjusted
    call(['cp', seqfile, 'rosetta'])
    os.chdir('rosetta')
    call([rosetta_make_fragments, seqfile, '-nosam', '-noporter'])
    os.chdir('..')
