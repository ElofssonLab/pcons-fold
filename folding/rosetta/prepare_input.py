import sys
import os
from subprocess import call

from reformat_contacts import reformat

rosettadir = "/home/x_mirmi/glob/rosetta"

seqfile = sys.argv[1]
cfile = sys.argv[2]
factor = sys.argv[3]

call(['mkdir', 'rosetta'])

### reformat the contact maps into constraint files used by rosetta
### creates: "cfile-<factor>.constraints"
reformat(seqfile, cfile, factor)
ros_cfile = cfile[:-4] + '-' + str(factor) + '.constraints'
call(['mv', ros_cfile, 'rosetta/'])

### run make_fragments.pl to create fragment database
### NOTICE: paths in make_fragments.pl need to be adjusted
call(['cp', seqfile, 'rosetta/'])
call(['cp', seqfile + '.ss2', 'rosetta/'])
os.chdir('rosetta')
call([rosettadir + '/rosetta_fragments/make_fragments.pl', seqfile, '-nosam', '-noporter', '-psipredfile', seqfile + '.ss2', '-id', '1xxxx'])


