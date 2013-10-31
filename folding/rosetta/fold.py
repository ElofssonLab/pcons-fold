import sys
import os
import math
from subprocess import call
from subprocess import Popen

### NOTE: Rosetta itself does not run in parallel, it uses a single core per run.
### Since it is non-deterministic (based on markov chain monte carlo) one can just
### run multiple instances at the same time and distribute them onto multiple cores.

#root = '/home/x_mirmi/pcons-fold/folding/rosetta/test/rosetta'
rosettadir = '/home/x_mirmi/glob/rosetta'

### default value for one triolith node
cores = 16

### change it if other value given
if sys.argv[1] == '-c':
    try:
        cores = int(sys.argv[2])
    except:
        print 'Number of cores -c must be an integer number, %r is not. Default is %s.' % (sys.argv[2], cores)
        sys.exit(1)
    del sys.argv[1:3]

root = os.path.abspath(sys.argv[1])
seqfile = os.path.abspath(sys.argv[2])
ros_cfile = os.path.abspath(sys.argv[3])

### The number of decoy structures produced in a single run depends on the 
### of cores. In total it should be at least 2000 to get native-like structures.
### (see PconsFold paper)
decoys_per_core = int(math.ceil(2000.0 / cores))


### 1..16 folders containing run 1..16
for core in range(1, cores + 1):
    os.chdir(root)
    optionstr = ''
    optionstr += '-in:file:fasta ' + seqfile + '\n'
    optionstr += '-in:file:frag3 ' + root + '/' + 'aa1xxxx03_05.200_v1_3\n'
    optionstr += '-in:file:frag9 ' + root + '/' + 'aa1xxxx09_05.200_v1_3\n'
    optionstr += '-use_filters true\n'
    optionstr += '-constraints:cst_file ' + ros_cfile + '\n'
    optionstr += '-database ' + rosettadir + '/rosetta_database/\n'
    optionstr += '-psipred_ss2 ' + seqfile + '.ss2/\n'
    optionstr += '-out:nstruct ' + str(decoys_per_core) + '\n'

    call(['mkdir', 'run_' + str(core)])
    os.chdir('/run_' + str(core))

    optionfile = open('abinitio.options', 'w')
    optionfile.write(optionstr)
    optionfile.close()

    Popen([rosettadir + '/rosetta_source/bin/AbinitioRelax.static.linuxgccrelease', '@abinitio.options'])
