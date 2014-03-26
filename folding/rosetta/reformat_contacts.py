import sys


def reformat(seqfile_name, infile_name, factor, outfile_name='', score=15):

    min_dist = 5
    factor = float(factor)

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
        if abs(res2 - res1) >= min_dist: 
            atm1 = 'CB'
            if seq[res1 - 1] == 'G':
                atm1 = 'CA'
            atm2 = 'CB'
            if seq[res2 - 1] == 'G':
                atm2 = 'CA'
            rosetta_lines.append('AtomPair %s %d %s %d FADE -10 19 10 %.2f 0' % (atm1, res1, atm2, res2, round((float(score) * -1.0) , 2)))
            count += 1
        if count > (seq_len * factor):
            break

    # write rosetta readable constraint file
    if not outfile_name:
        outfile_name = infile_name + '-' + str(factor) + '.constraints'
    outfile = open(outfile_name, 'w')
    for line in rosetta_lines:
        outfile.write('%s\n' % line)



if __name__ == "__main__":

    seqfile_name = sys.argv[1]
    infile_name = sys.argv[2]

    # use L * factor highest scoring constraints
    factor = float(sys.argv[3])
    score = sys.argv[4]

    reformat(seqfile_name, infile_name, factor, score)
