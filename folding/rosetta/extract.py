#!/usr/bin/env python
import sys
import os
from subprocess import call

import evaluate_scores as ev

from localconfig import *


def main(seqfile, n_cores=1, n_models=10, relax_flag=True, rundir_postfix='', nametag='', chain='A'):

    rundir = os.path.dirname(os.path.abspath(seqfile)) + '/'

    ### determine IDs of top scoring decoys
    all_scores_sorted = ev.get_best_of_all_runs(n_models, n_cores, rundir + rundir_postfix)

    ### extract given IDs
    ev.extract_structures(all_scores_sorted)

    ### and relax the models if given (RECOMMENDED)
    if relax_flag:
        ev.relax_structures(all_scores_sorted)

    ### calculate similarity scores to native structure, if given
    ### and store them in "TMscores.txt"
    native_fname = rundir + 'native.pdb'
    if os.path.exists(native_fname):
        scorestr = ev.compare_to_native(all_scores_sorted, relax_flag, False, native_fname=native_fname, name=nametag, chain=chain)
        scorefile = open(rundir + rundir_postfix + '/TMscores.txt', 'w')
        scorefile.write(scorestr)
        scorefile.close()
     

if __name__ == "__main__":

    ### default value = relax top-ranked structures
    relax_flag = True
    ### default: extract top 10 ranked models
    n_models = 10
    
    rundir_postfix = ''
    nametag = ''
    chain = 'A'

    ### parse parameters
    if '-c' in sys.argv:
        idx = sys.argv.index('-c')
        try:
            n_cores = int(sys.argv[idx+1])
        except:
            print 'Number of cores -c must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_cores)
            sys.exit(1)
        del sys.argv[idx]
        del sys.argv[idx]

    if '-m' in sys.argv:
        idx = sys.argv.index('-m')
        try:
            n_models = int(sys.argv[idx+1])
        except:
            print 'Number of models to extract -m must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_models)
            sys.exit(1)
        del sys.argv[idx]
        del sys.argv[idx]

    if '--norelax' in sys.argv:
        idx = sys.argv.index('--norelax')
        relax_flag = False
        del sys.argv[idx]

    if '--rundir_postfix' in sys.argv:
        idx = sys.argv.index('--rundir_postfix')
        rundir_postfix = sys.argv[idx+1]
        del sys.argv[idx]
        del sys.argv[idx]

    if '--nametag' in sys.argv:
        idx = sys.argv.index('--nametag')
        nametag = sys.argv[idx+1]
        del sys.argv[idx]
        del sys.argv[idx]

    if '--chain' in sys.argv:
        idx = sys.argv.index('--chain')
        chain = sys.argv[idx+1]
        del sys.argv[idx]
        del sys.argv[idx]

    seqfile = os.path.abspath(sys.argv[1])

    main(seqfile, n_cores=n_cores, n_models=n_models, relax_flag=relax_flag, rundir_postfix=rundir_postfix, nametag=nametag, chain=chain)
