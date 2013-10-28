#!/usr/bin/env python
import sys, os, subprocess

if __name__ == '__main__':
	print 'Please do not run me! Use predict.py or predictAll.py'
	print '\n\tYours sincerely,\n\n\t', sys.argv[0]
	sys.exit(0)

sys.stderr.write("""
****************************************************************************
          PconsC : Combination of direct information methods and alignments 
                   improves contact prediction
****************************************************************************

If you use PconsC for contact prediction, please cite the following papers:

Skwark M.J., Abdel-Rehim A. and Elofsson A.  
"PconsC : Combination of direct information methods and alignments improves 
          contact prediction"

Ekeberg M, Lovkvist C, Lan Y., Weigt M. and Aurell E.
"Improved contact prediction in proteins: Using pseudolikelihoods to infer 
          Potts models"
doi: 10.1103/PhysRevE.87.012707

Jones D.T, Buchan D.W.A., Cozzetto D., Pontii M. 
"PSICOV: precise structural contact prediction using sparse inverse covariance 
         estimation on large multiple sequence alignments"
doi: 10.1093/bioinformatics/btr638

Johnson L., Eddy S. and Portugaly E.
"Hidden Markov model speed heuristics and iterative HMM search procedure"
doi:10.1186/1471-2105-11-431

Remmert M., Biegert A., Hauser A. and Soding J.     
"HHblits: lightning-fast iterative protein sequence searching by HMM-HMM
          alignment"
doi:10.1038/nmeth.1818 
-----------------------------------------------------------------------------

""")

# Directory where the distributable package is located
# e.g. root = '/afs/pdc.kth.se/home/m/mjs/volume4/correlated-bench/'
root = os.path.dirname(os.path.abspath(sys.argv[0])) + '/'

# Maximum amount of cores to use
cores = 4

# Enable work-around for PSICOV not handling low complexity alignments?
# Warning: enable ONLY if you are certain that it is what you want to do!
# To enable, change False into True
psicovfail = False

# Path to formatdb formatted sequence database (e.g. Uniref90, nr90 etc.)
# We recommend UniRef100
# e.g. jackhmmerdb = '/home/mjs/db/uniref/uniref100.fasta'
jackhmmerdb = ''


# Path to HHblits database
# e.g. hhblitsdb = '/home/mjs/db/hhpred/new/nr20_12Aug11'
hhblitsdb = ''

# Path to MATLAB executable
# e.g. matlab = '/afs/pdc.kth.se/pdc/vol/matlab/r2012a/bin/matlab'
matlab = ''

# Path to executable files
jackhmmer = ''
hhblits = ''
psicov = ''

# These are included. Should not need changing.
scriptpath = root + 'scripts'
trim = root + 'scripts/trim.py'
trim2 = root + 'scripts/trimToFasta.py'

# Reformat script scavenged from HHsuite. Please cite the HHblits paper!
reformat = root + 'scripts/reformat.pl'

# Download plmDCA from http://plmdca.csc.kth.se/ and put it into scripts/plmDCA_symmetric_v2 directory
# e.g. plmdca = root + 'scripts/plmDCA_symmetric_v2/plmDCA_symmetric.m'

plmdca = root + 'scripts/plmDCA_symmetric_v2/plmDCA_symmetric.m'
if not os.path.exists(plmdca):
	sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
	sys.stderr.write('Download plmDCA from http://plmdca.csc.kth.se/ and put it into scripts/plmDCA_symmetric_v2 directory!\n')
	sys.exit(1)

if len(jackhmmerdb) < 2:
	sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
	sys.stderr.write('jackhmmer database NOT SET!\n')
	sys.exit(1)

if not os.path.exists(jackhmmerdb):
	sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
	sys.stderr.write('jackhmmer database (' + jackhmmerdb + ') DOES NOT EXIST!\n')
	sys.exit(1)

if len(hhblitsdb) < 2:
	sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
	sys.stderr.write('HHblits database NOT SET!\n')
	sys.exit(1)

if not os.path.exists(hhblitsdb+'_a3m_db'):
	sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
	sys.stderr.write('HHblits database (' + hhblitsdb + ') DOES NOT EXIST!\n')
	sys.exit(1)

try:
	f = open(os.devnull, "w") 
	x  = subprocess.call([matlab, '-h'], stdout=f, stderr=f)
	f.close()
except:
	sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
	sys.stderr.write('Chosen MATLAB binary does not seem to work!\n')
	sys.exit(1)

try:
	f = open(os.devnull, "w") 
	x  = subprocess.call([jackhmmer, '-h'], stdout=f, stderr=f)
	f.close()
except Exception as e:
	sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
	sys.stderr.write('Chosen jackhmmer binary does not seem to work!\n')
	sys.exit(1)

try:
	f = open(os.devnull, "w") 
	x  = subprocess.call([hhblits, '-h'], stderr=f, stdout=f)
	f.close()
	pass
except:
	sys.stderr.write('*****************\n   ERROR!\n*****************\n\n')
	sys.stderr.write('Chosen HHblits binary does not seem to work!\n')
	sys.exit(1)


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
