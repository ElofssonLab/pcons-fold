import sys
import os
from subprocess import call

from reformat_contacts import reformat

rosettadir = "home/x_mirmi/glob/rosetta"

seqfile = sys.argv[1]
cfile = sys.argv[2]
factor = sys.argv[3]

### reformat the contact maps into constraint files used by rosetta
### creates: "cfile-<factor>.constraints"
reformat(seqfile, cfile, factor)

### run make_fragments.pl to create fragment database
### NOTICE: paths in make_fragments.pl need to be adjusted
call(['mkdir', 'rosetta'])
call(['mv', seqfile + '-' + str(factor) + '.constraints', 'rosetta/'])
os.chdir('rosetta')
call([rosettadir + '/rosetta_fragments/make_fragments.pl', seqfile, '-nosam', '-noporter', '-psipredfile', seqfile + '.ss', '-id', '1xxxx'])


