#!/usr/bin/env python
import os
import sys
import multiprocessing
import subprocess

if __name__ == '__main__':
    print 'Please do not run me! Use run_pipeline.py'
    print '\n\tYours sincerely,\n\n\t', sys.argv[0]
    sys.exit(0)

sys.stderr.write("""
****************************************************************************
     PconsFold : Improved contact predictions improve protein models 
****************************************************************************

If you use PconsFold for protein structure prediction, 
please cite the following publication:

-----------------------------------------------------------------------------

""")

# Directory where the distributable package is located
root = os.path.dirname(os.path.abspath(sys.argv[0])) + '/'

# Look if PconsC or Rosetta dependencies need to be checked
rosetta_flag = False
if 'rosetta' in root:
    rosetta_flag = True
    root = '/'.join(root.split('/')[:-3]) + '/pconsc/'



########################################################
### Please adjust the following paths to your system ###
########################################################

### Path to root folder of your Rosetta installation ###
# REQUIRES: Rosetta 3.5 or weekly build
rosettadir = '/home/mirco_local/scratch/apps/rosetta_2013wk42_bundle'

### Jackhmmer executable ###
jackhmmer = root + 'dependencies/hmmer-3.0/src/jackhmmer'

### HHblits executable ###
hhblits = root + 'dependencies/hhsuite-2.0.16/bin/hhblits'

### PSICOV executable ###
psicov = root + 'dependencies/psicov-1.11/psicov'

### MATLAB executable ###
# Please set this variable to None if you don't have access to matlab. 
# PconsFold will then try to use the compiled version. 
matlab = '/sw/apps/matlab/x86_64/8.1/bin/matlab'
#matlab = None

### Path to MATLAB compiler ###
# Only needed if matlab is not available.
matlabdir = '/software/apps/mcr/2012b/build01/v80/' 

### Path to TM-score ###
# only needed if result should be compared to native structure
tmscore_binary = '/home/x_mirmi/pcons-fold/folding/rosetta/dependencies/TMscore/run_TMscore.sh'



########################################################
###  Please do not change anything below this line   ###
########################################################


# Rosetta paths
rosetta_db_dir = rosettadir + '/main/database'
rosetta_binary_dir = rosettadir + '/main/source/bin'
rosetta_make_fragments = rosettadir + '/tools/fragment_tools/make_fragments.pl'
rosetta_abinitiorelax = rosetta_binary_dir + '/AbinitioRelax.linuxgccrelease'
rosetta_extract = rosetta_binary_dir + '/extract_pdbs.linuxgccrelease'
rosetta_relax = rosetta_binary_dir + '/relax.linuxgccrelease'


# Paths to included scripts
trim2jones = root + 'scripts/a3mToJones.py'
trim2trimmed = root + 'scripts/a3mToTrimmed.py'
#trim = root + 'scripts/trim.py'
#trim2 = root + 'scripts/trimToFasta.py'

# Reformat script scavenged from HHsuite. Please cite the HHblits paper!
reformat = root + 'scripts/reformat.pl'

# Maximum amount of cores to use per default
n_cores = multiprocessing.cpu_count()

# Enable work-around for PSICOV not handling low complexity alignments
psicovfail = True

# Adjust plmdca path to either standalone or compiled, 
# depending on presence of matlab.
if matlab:
    plmdca = None # matlab licence present: do not use compiled version
    plmdcapath = root + 'dependencies/plmDCA_symmetric_v2'
else:
    plmdca = root + 'dependencies/plmdca_standalone/2012/build01/bin/plmdca'
    plmdcapath = None



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
            sys.exit(1)
    else:
        sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
        sys.stderr.write('You must set one of plmdca or matlab in localconfig.py!\n')
        sys.exit(1)


sys.stderr.write('\nDependencies OK.\n')
