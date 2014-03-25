import sys
import os
import evaluate_scores as ev

from localconfig import *

### default value = detected number of cores
cores = 16
### default value = relax top-ranked structures
relax_flag = True

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

if '--norelax' in sys.argv:
    idx = sys.argv.index('--norelax')
    relax_flag = False
    del sys.argv[idx]

num = int(sys.argv[1])

### determine IDs of top scoring decoys
all_scores_sorted = ev.get_best_of_all_runs(num, cores, rundir)

### extract given IDs
ev.extract_structures(all_scores_sorted)

### and relax the models if given (RECOMMENDED)
if relax_flag:
    ev.relax_structures(all_scores_sorted)

### calculate similarity scores to native structure, if given
### and store them in "TMscores.txt"
if os.path.exists(rundir + '/native.pdb'):
    scorestr = ev.compare_to_native(all_scores_sorted, relax_flag, False)
    scorefile = open(rundir + '/TMscores.txt', 'w')
    scorefile.write(scorestr)
    scorefile.close()
 



