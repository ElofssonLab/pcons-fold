#!/usr/bin/python

import pickle, sys, os, random

count = 0
forest = pickle.load(open(os.path.dirname(os.path.abspath(sys.argv[0])) + '/layer0.forest'))

if len(sys.argv) != 20:
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
	print ' NetSurf RSA'
	print ' PSIPRED SS'
	print ' Output file', len(sys.argv)
	sys.exit(1)

def parseNetSurfP(f):
        netSurfdict = {}
        for l in open(f).readlines():
                al = []
                x = l.split()
                if l.find('#') == 0:
                        continue
#               for y in [4,6,7, 8, 9]:
#                       al.append(float(x[y]) )
                if x[0] == 'B':
                        al = [0]
                else:
                        al = [1]
                netSurfdict[ int(x[3] )] = al
        return netSurfdict


def parsePSIPREDhoriz(hfile):
        SS = ''
        conf = ''
        SSdict = {}
        for line in hfile:
                line_arr = line.strip().split(' ')
                if line_arr[0] == 'Pred:' and len(line_arr) > 1:
                        SS += line_arr[1]
                if line_arr[0] == 'Conf:' and len(line_arr) > 1:
                        conf += line_arr[1]
        for i in range(len(conf)):
                if SS[i] == 'H':
                        SSdict[i+1] = [int(conf[i]), 0,0]
                elif SS[i] == 'E':
                        SSdict[i+1] = [0, int(conf[i]), 0]
                else:
                        SSdict[i+1] = [0,0,int(conf[i])]
        return SSdict


def parsePSIPRED(f):
        x = open(f).read().split('\n')
        conf = x[3]
        SS = x[1]
        SSdict = {}
        for i in range(len(conf)):
                if SS[i] == 'H':
                        SSdict[i+1] = [int(conf[i]), 0,0]
                elif SS[i] == 'E':
                        SSdict[i+1] = [0, int(conf[i]), 0]
                else:
                        SSdict[i+1] = [0,0,int(conf[i])]
        return SSdict


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

selected = set()
contacts = {}
X = []
Y = []
maxres = -1
acceptable = []
accessibility = parseNetSurfP(files[16])
#SSdict = parsePSIPRED(files[17])
SSdict = parsePSIPREDhoriz(files[17])

for index in range(16):
	contacts[index] = {}
	d = files[index]
	r = []
	if not os.path.exists(d):
		sys.stderr.write(d + ' does not exist!\n')
		sys.exit(1)
	infile = open(d).readlines()
	if len(infile) < 1 and index not in (6, 14):
		sys.stderr.write(d + ' is empty! Performance MAY be affected\n')
	else:
		acceptable.append(index)
	for m in infile:
		if d.find('psicov') > -1:
			x = m.split()
			if len(x) != 5:
				print d + ' has wrong format!'
				sys.exit(1)
			c = 4
		elif d.find('plmdca') > -1:
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
		if abs(aa1 - aa2) < 6:
			continue
		if x[c].find('nan') > -1:
			score = -3
		else:
			score = float(x[c])
		contacts[index][(aa1, aa2)] = score

clist = []
for c in contacts[0].keys():
	q = [ c ]
	for i in contacts.keys():
		try:
			q.append( contacts[i][c] )
		except:
			q.append( -3 )
	clist.append(q)


for j in clist:
	selected.add(j[0])


outf = open(files[18], 'w')
for s in sorted(list(selected)):
	q = []
	for index in range(16):
		try:
			q.append(contacts[index][s])
		except:
			q.append(-3)
	
	for i in range(-4, 5):
		try:
			q.extend(SSdict[s[0]] )
		except:
			q.extend([0,0,0])
	
	for i in range(-4, 5):
		try:
			q.extend(SSdict[s[1]] )
		except:
			q.extend([0,0,0])

	
	for i in range(-4, 5):
		try:
			q.extend(accessibility[s[0]] )
		except:
#			q.extend([0,0,0,0,0])
			q.extend([-1])
	
	for i in range(-4, 5):
		try:
			q.extend(accessibility[s[1]] )
		except:
#			q.extend([0,0,0,0,0])
			q.extend([-1])

	outf.write('%d %d %6.4f\n' % (s[0], s[1], predict(q, forest) ) )

outf.close()

