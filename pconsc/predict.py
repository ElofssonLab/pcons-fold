#!/usr/bin/env python

"""Script to predict contacts, given the 16 input files"""

import pickle, sys
import os

if len(sys.argv) != 17:
	print 'Usage: ' + sys.argv[0] + ' <output files>'
	print 'Output files need to come in *order*!'
	print 'That is:'
	print ' JackHMMER 1e-4 Psicov'
	print ' JackHMMER 1e-4 plmDCA'
	print ' JackHMMER 1e-0 Psicov'
	print ' JackHMMER 1e-0 plmDCA'
	print ' JackHMMER 1e-10 Psicov'
	print ' JackHMMER 1e-10 plmDCA'
	print ' JackHMMER 1e-40 Psicov'
	print ' JackHMMER 1e-40 plmDCA'
	print ' HHblits 1e-4 Psicov'
	print ' HHblits 1e-4 plmDCA'
	print ' HHblits 1e-0 Psicov'
	print ' HHblits 1e-0 plmDCA'
	print ' HHblits 1e-10 Psicov'
	print ' HHblits 1e-10 plmDCA'
	print ' HHblits 1e-40 Psicov'
	print ' HHblits 1e-40 plmDCA'
	sys.exit(1)


def predict(X, forest):
	probability = []
	for t in range(len(forest)):
		tree = forest[t]
		while len(tree) > 2:
			if X[tree[0][0]] <= tree[0][1]:
				tree = tree[1]
			else:
				tree = tree[2]
		probability.append(tree[1]/float(tree[0] + tree[1]))
	return sum(probability)/len(probability)

files = sys.argv[1:]
forest = pickle.load(open(os.path.dirname(os.path.abspath(sys.argv[0])) + '/forest.dat'))

selected = set()
contacts = {}
X = []
Y = []
maxres = -1
for index in range(16):
	contacts[index] = {}
	d = files[index]
	r = []
	if not os.path.exists(d):
		sys.stderr.write(d + ' does not exist!\n')
		continue
#		sys.exit(1)
	infile = open(d).readlines()
	if len(infile) < 1:
		sys.stderr.write(d + ' is empty!\n')
		continue
#		sys.exit(1)
	for m in infile:
		if index % 2 == 0:
			x = m.split()
			if len(x) != 5:
				print d + ' has wrong format!'
				sys.exit(1)
			c = 4
		else:
			x = m.split(',')
			if len(x) != 3:
				print d + ' has wrong format!'
				sys.exit(1)
			c = 2
		aa1 = int(x[0])
		aa2 = int(x[1])
		if aa1 > maxres:
			maxres = aa1
		if aa2 > maxres:
			maxres = aa2	
		if abs(aa1 - aa2) < 5:
			continue
		if x[c].find('nan') > -1:
			score = 0
		else:
			score = float(x[c])
		selected.add( (aa1, aa2) )
		contacts[index][(aa1, aa2)] = score

for s in sorted(list(selected)):
	q = []
	for index in range(16):
		try:
			q.append(contacts[index][s])
		except:
			q.append(0)

	if len(q) == 16:
		X.append(q)
		Y.append(s)

for l in range(len(Y)):
	(aa1, aa2) = (Y[l][0], Y[l][1])
	print '%d %d %6.4f' % (aa1, aa2, predict(X[l], forest))


