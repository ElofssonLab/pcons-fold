import sys
import os
from subprocess import call

import reformat_contacts
from localconfig import *

nohoms_flag = False

if '--nohoms' in sys.argv:
    idx = sys.argv.index('--nohoms')
    nohoms_flag = False
    del sys.argv[idx]

seqfile = sys.argv[1]
cfile = sys.argv[2]
factor = sys.argv[3]

saved_path = os.getcwd()
rundir = os.path.dirname(os.path.abspath(seqfile)) + '/'

if not os.path.exists(rundir + 'rosetta'):
    call(['mkdir', rundir + 'rosetta'])

### reformat the contact maps into constraint files used by rosetta
### creates: "cfile-<factor>.constraints"
ros_cfile = cfile + '-' + str(factor) + '.constraints'

reformat_contacts.reformat(seqfile, cfile, factor, outfile_name=ros_cfile)
call(['mv', ros_cfile, rundir + 'rosetta'])

### run make_fragments.pl to create fragment database
### NOTICE: paths in make_fragments.pl need to be adjusted
call(['cp', seqfile, rundir + 'rosetta'])
os.chdir(rundir + 'rosetta')
if nohoms_flag:
    call([rosetta_make_fragments, os.path.basename(seqfile), '-nosam', '-noporter', '-nohoms'])
else:
    call([rosetta_make_fragments, os.path.basename(seqfile), '-nosam', '-noporter'])

os.chdir(saved_path)
