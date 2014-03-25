import sys
import os
import math
from subprocess import call
from subprocess import Popen

from localconfig import *

### default value
n_cores = 16
n_decoys = 2000

### parse parameters
if '-c' in sys.argv:
    idx = sys.argv.index('-c')    
    try:
        n_cores = int(sys.argv[idx+1])
    except:
        print 'Number of cores -c must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_cores)
        sys.exit(1)
    del sys.argv[idx]
    del sys.argv[idx+1]

if '-n' in sys.argv:
    idx = sys.argv.index('-n')    
    try:
        n_decoys = int(sys.argv[idx+1])
    except:
        print 'Number of decoys -n must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_decoys)
        sys.exit(1)
    del sys.argv[idx]
    del sys.argv[idx+1]

seqfile = os.path.abspath(sys.argv[1])
ros_cfile = os.path.abspath(sys.argv[2])

### The number of decoy structures produced in a single run depends on the 
### of cores. In total it should be at least 2000 to get native-like structures.
### (see PconsFold paper)
decoys_per_core = int(math.ceil(float(n_decoys) / n_cores))


### n_cores folders containing run 1..n_cores
plist = []
for core in range(1, n_cores + 1):
    optionstr = ''
    optionstr += '-in:file:fasta ' + seqfile + '\n'
    optionstr += '-in:file:frag3 ' + rundir + '/aa1xxxx03_05.200_v1_3\n'
    optionstr += '-in:file:frag9 ' + rundir + '/aa1xxxx09_05.200_v1_3\n'
    optionstr += '-use_filters true\n'
    optionstr += '-constraints:cst_file ' + ros_cfile + '\n'
    optionstr += '-database ' + rosettadir + '/rosetta_database/\n'
    optionstr += '-psipred_ss2 ' + seqfile + '.ss2/\n'
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
