#!/usr/bin/env python
import sys, os
import shutil
import subprocess
import multiprocessing

from localconfig import *
from pconsc import predict_all
from folding.rosetta import prepare_input
from folding.rosetta import fold
from folding.rosetta import extract

sys.stderr.write("""
****************************************************************************
     PconsFold : Improved contact predictions improve protein models 
****************************************************************************

If you use PconsFold for protein structure prediction, 
please cite the following publication:

-----------------------------------------------------------------------------

""")

### parse parameters
if len(sys.argv) < 4:
    sys.stderr.write('Usage: ./%s [-c n_cores] [-n n_decoys] [-m n_models]\n' % sys.argv[0].strip('./'))
    sys.stderr.write('            [-f factor] [--norelax] [--nohoms]\n')
    sys.stderr.write('            <hhblits db> <jackhmmer db> <sequence file>\n')
    sys.exit(0)

sys.stderr.write('\nTesting dependencies...\n')

if rosetta_flag:

    ### Check Rosetta ###
    try:
        f = open(os.devnull, "w")
        x  = subprocess.call([rosetta_make_fragments, '-h'], stderr=f, stdout=f)
        f.close()
        pass
    except:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('There might be something wrong with your Rosetta installation in:\n')
        sys.stderr.write(rosettadir + '\n')
        sys.stderr.write('Please check the path to the Rosetta root directory\n')
        sys.stderr.write('and use Rosetta 3.5 or higher (weekly).\n')
        sys.stderr.write('Please ensure that the following Rosetta executable\n')
        sys.stderr.write('is present and working:\n')
        sys.stderr.write(rosetta_make_fragments + '\n')
        sys.exit(1)
    try:
        f = open(os.devnull, "w")
        x  = subprocess.call([rosetta_abinitiorelax, '-h'], stderr=f, stdout=f)
        f.close()
        pass
    except:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('There might be something wrong with your Rosetta installation in:\n')
        sys.stderr.write(rosettadir + '\n')
        sys.stderr.write('Please check the path to the Rosetta root directory\n')
        sys.stderr.write('and use Rosetta 3.5 or higher (weekly).\n')
        sys.stderr.write('Please ensure that the following Rosetta executable\n')
        sys.stderr.write('is present and working:\n')
        sys.stderr.write(rosetta_abinitiorelax + '\n')
        sys.exit(1)
    try:
        f = open(os.devnull, "w")
        x  = subprocess.call([rosetta_extract, '-h'], stderr=f, stdout=f)
        f.close()
        pass
    except:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('There might be something wrong with your Rosetta installation in:\n')
        sys.stderr.write(rosettadir + '\n')
        sys.stderr.write('Please check the path to the Rosetta root directory\n')
        sys.stderr.write('and use Rosetta 3.5 or higher (weekly).\n')
        sys.stderr.write('Please ensure that the following Rosetta executable\n')
        sys.stderr.write('is present and working:\n')
        sys.stderr.write(rosetta_extract + '\n')
        sys.exit(1)
    try:
        f = open(os.devnull, "w")
        x  = subprocess.call([rosetta_relax, '-h'], stderr=f, stdout=f)
        f.close()
        pass
    except:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('There might be something wrong with your Rosetta installation in:\n')
        sys.stderr.write(rosettadir + '\n')
        sys.stderr.write('Please check the path to the Rosetta root directory\n')
        sys.stderr.write('and use Rosetta 3.5 or higher (weekly).\n')
        sys.stderr.write('Please ensure that the following Rosetta executable\n')
        sys.stderr.write('is present and working:\n')
        sys.stderr.write(rosetta_relax + '\n')
        sys.exit(1)
else:

    ### Check Jackhmmer ###
    try:
        f = open(os.devnull, "w")
        x  = subprocess.call([jackhmmer, '-h'], stdout=f, stderr=f)
        f.close()
    except Exception as e:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('Chosen jackhmmer binary does not seem to work!\n')
        sys.exit(1)

    ### Check HHblits ###
    try:
        f = open(os.devnull, "w")
        x  = subprocess.call([hhblits, '-h'], stderr=f, stdout=f)
        f.close()
        pass
    except:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('Chosen HHblits binary does not seem to work!\n')
        sys.exit(1)

    ### Check PSICOV ###
    try:
        f = open(os.devnull, "w")
        x  = subprocess.call([psicov, root + '/extras/psicovtest.fas'], stdout=f, stderr=f)
        f.close()
    except Exception as e:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('Chosen PSICOV binary does not seem to work!\n')
        sys.exit(1)

    if x == 255 and not psicovfail:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('Your version of PSICOV refuses to handle low-complexity alignments.\n')
        sys.stderr.write('We recommend patching the PSICOV code to allow this. See 00README\n')
        sys.stderr.write('If you _really_ do not want to do that, please change psicovfail flag in \n')
        sys.stderr.write(os.path.abspath(sys.argv[0]) + ' to True.\n')
        sys.stderr.write('This will (most probably) affect the prediction performance.\n')
        sys.exit(1)


    ### Check plmDCA ###
    if plmdca:
        try:
            f = open(os.devnull, "w")
            x  = subprocess.call([plmdca, '-h'], stdout=f, stderr=f)
            f.close()
        except Exception as e:
            sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
            sys.stderr.write('Chosen plmdca binary does not seem to work!\n')
            sys.exit(1)
    elif matlab:
        try:
            f = open(os.devnull, "w")
            x  = subprocess.call([matlab, '-h'], stdout=f, stderr=f)
            f.close()
        except:
            sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
            sys.stderr.write('Chosen MATLAB binary does not seem to work!\n')
            sys.stderr.write('You can get MCR \n')
            sys.stderr.write('http://www.mathworks.se/products/compiler/mcr/\n')
            sys.exit(1)
    else:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('You must set one of plmdca or matlab in localconfig.py!\n')
        sys.exit(1)


sys.stderr.write('Dependencies OK.\n')


nohoms_flag = False
relax_flag = True
factor = 1.0
n_models = 10
n_decoys = 2000

if '-c' in sys.argv:
    idx = sys.argv.index('-c')
    try:
        n_cores = int(sys.argv[idx+1])
    except:
        print 'Number of cores -c must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_cores)
        sys.exit(1)
    del sys.argv[idx]
    del sys.argv[idx]

if '-n' in sys.argv:
    idx = sys.argv.index('-n')
    try:
        n_decoys = int(sys.argv[idx+1])
    except:
        print 'Number of decoys -n must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], n_decoys)
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

if '-f' in sys.argv:
    idx = sys.argv.index('-f')
    try:
        factor = int(sys.argv[idx+1])
    except:
        print 'Factor of sequence length (determining number of constraints to be used during folding) -f must be float, %r is not. Default is %s.' % (sys.argv[idx+1], factor)
        sys.exit(1)
    del sys.argv[idx]
    del sys.argv[idx]

if '--nohoms' in sys.argv:
    idx = sys.argv.index('--nohoms')
    nohoms_flag = True
    del sys.argv[idx]

if '--norelax' in sys.argv:
    idx = sys.argv.index('--norelax')
    relax_flag = False
    del sys.argv[idx]


hhblitsdb = os.path.abspath(sys.argv[1])
jackhmmerdb = os.path.abspath(sys.argv[2])
seqfile = os.path.abspath(sys.argv[3])
contactfile = seqfile + '.pconsc.out'
constraintfile = contactfile + '-' + str(factor) + '.constraints'

shutil.copyfile(root + '../localconfig.py', root + 'localconfig.py')
shutil.copyfile(root + '../localconfig.py', root + '../folding/rosetta/localconfig.py')

predict_all.main(hhblitsdb, jackhmmerdb, seqfile, n_cores=n_cores)

rundir_postfix = 'rosetta'
prepare_input.main(seqfile, contactfile, factor=factor, nohoms_flag=nohoms_flag)
fold.main(seqfile, constraintfile, n_cores=n_cores, n_decoys=n_decoys, rundir_postfix=rundir_postfix)    
extract.main(seqfile, n_cores=n_cores, n_models=n_models, relax_flag=relax_flag, rundir_postfix=rundir_postfix)


### collect the results in seperate folder
call(['mkdir', rundir + 'rosetta_results'])
call('mv %s/%s/*.run_*.*.pdb %s/rosetta_results' % (rundir, rundir_postfix, rundir), shell=True)

if os.path.exists(rundir + 'native.pdb'):
    call(['mv', rundir + rundir_postfix + '/TMscores.txt', rundir + 'rosetta_results'])

