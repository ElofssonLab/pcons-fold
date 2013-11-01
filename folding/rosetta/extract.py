import sys
import os
import evaluate_scores as ev

from localconfig import *

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

num = int(sys.argv[1])
relax_flag = bool(int(sys.argv[2]))

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
 



