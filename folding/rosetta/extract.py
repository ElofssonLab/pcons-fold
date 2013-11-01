import sys
import os
import evaluate_scores as ev

num = int(sys.argv[1])
cores = int(sys.argv[2])
rundir = sys.argv[3]
relax_flag = bool(int(sys.argv[4]))

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
 



