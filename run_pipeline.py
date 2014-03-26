import sys, os
import shutil
from subprocess import call

### parse parameters
if len(sys.argv) < 4:
    print 'Usage:', sys.argv[0], '[-c n_cores] [-n n_decoys] <hhblits db> <jackhmmer db> <sequence file>'
    sys.exit(0)

if '-c' in sys.argv:
    idx = sys.argv.index('-c')
    try:
        n_cores = int(sys.argv[idx+1])
    except:
        print 'Number of cores -c must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_cores)
        sys.exit(1)
    del sys.argv[idx +1]
    del sys.argv[idx]

if '-n' in sys.argv:
    idx = sys.argv.index('-n')
    try:
        n_decoys = int(sys.argv[idx+1])
    except:
        print 'Number of decoys -n must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_decoys)
        sys.exit(1)
    del sys.argv[idx +1]
    del sys.argv[idx]

hhblitsdb = sys.argv[1]
jackhmmerdb = sys.argv[2]
seqfile = sys.argv[3]


def check_output(command, ok_to_fail=False):
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if ok_to_fail or not p.returncode:
        return stdout
    msg = "ERROR: Command %s exited with nonzero returncode %s. Stdout available in %s, stderr in %s"
    err = seqfile + '.fail.stderr'
    out = seqfile + '.fail.stdout'
    print >> open(err, 'w'), stderr
    print >> open(out, 'w'), stdout
    print >> sys.stderr, msg % (command, p.returncode, out, err)
    sys.exit(1)



shutil.copyfile('localconfig.py', 'pconsc/localconfig.py')
shutil.copyfile('localconfig.py', 'folding/rosetta/localconfig.py')

python ../pcons-fold/pconsc/predictAll_1.0.py ../databases/hhsuite_db/uniprot20/uniprot20_2012_10_klust20_dc_2012_12_10 ../databases/uniref/current_release/uniref90.fasta  BPT1_BOVIN.fa

python ~/scratch/pcons-fold/folding/rosetta/prepare_input.py BPT1_BOVIN.fa BPT1_BOVIN.fa.pconsc.out 1.0

python ../pcons-fold/folding/rosetta/fold.py -c 2 -n 4 BPT1_BOVIN.fa rosetta/BPT1_BOVIN.fa.pconsc.out-1.0.constraints

python ../pcons-fold/folding/rosetta/extract.py 4 -c 2
