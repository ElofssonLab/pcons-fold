import sys

# command line input
seqfile_name = sys.argv[1]
infile_name = sys.argv[2]
#min_dist = int(sys.argv[3])
min_dist = 5

# set length of a single repeat unit:
rep_len = 33

# use L * factor highest scoring constraints
factor = float(sys.argv[3])
score = sys.argv[4]

# scale value x from [min_x, max_x] to [0,1]
# NOT USED
def scale(x, min_x, max_x):
    range_x = (max_x - min_x)
    #print range_x
    return (x - min_x) / range_x


# read sequence needed for checking if res = glycine
# to assign CA/CB
seq = ''
seqfile = open(seqfile_name, 'r')
for line in seqfile:
    if line[0] != '>':
        seq += line.strip()
seqfile.close()
seq_len = len(seq)


# guessing separator of constraint file
test_line = open(infile_name,'r').readline()
if len(test_line.split(',')) != 1:
    sep = ','
elif len(test_line.split(' ')) != 1:
    sep = ' '
else:
    sep = '\t'


# read constraint file
infile = open(infile_name, 'r')
old_constraints = []
for line in infile:
    old_constraints.append(tuple(line.strip().split(sep)))
    

# reformat to rosetta constraints
old_constraints.sort(key=lambda x: float(x[-1]), reverse=True)
min_score = float(old_constraints[-1][-1])
max_score = float(old_constraints[0][-1])
rosetta_lines = []
count = 0
for constr in old_constraints:
    res1 = int(constr[0])
    res2 = int(constr[1])
    if abs(res2 - res1) >= min_dist: # and abs(res2 - res1) < rep_len * 1.5:
        atm1 = 'CB'
        if seq[res1 - 1] == 'G':
            atm1 = 'CA'
        atm2 = 'CB'
        if seq[res2 - 1] == 'G':
            atm2 = 'CA'
        #score = scale(float(constr[2]), min_score, max_score)
        #score = float(constr[-1])
        #score = score * -2.0
        #score = -1.0 * score
        #score = -15.0
        rosetta_lines.append('AtomPair %s %d %s %d FADE -10 19 10 %.2f' % (atm1, res1, atm2, res2, round((float(score) * -1.0) , 2)))
        #rosetta_lines.append('AtomPair %s %d %s %d BOUNDED 1.5 8 1 0.5 PREDICTED' % (atm1, res1, atm2, res2))
        count += 1
    if count > (seq_len * factor):
        break


# write rosetta readable constraint file
outfile_name = '.'.join(infile_name.split('.')[0:-1]) + '-' + str(factor) + '.constraints'
#outfile_name = '.'.join(infile_name.split('.')[0:-1]) + '-s' + score + '.constraints'
#outfile_name = '.'.join(infile_name.split('.')[0:-1]) + '-BOUNDED.constraints'
outfile = open(outfile_name, 'w')
for line in rosetta_lines:
    outfile.write('%s\n' % line)
