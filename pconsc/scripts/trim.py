#!/usr/bin/env python

import sys, os

infilef = sys.argv[1]

if len(sys.argv) > 2:
        lookfor = sys.argv[2]
else:
	lookfor = 'target'

onetonumber = 'ARNDCEQGHILKMFPSTWYV-'

infile = open(infilef).read().split('\n')

sequences = {}
key = 'n/a'
trimto = 'target'

for l in infile:
	if l.find('>') == 0:
		if key == trimto:
			break
		key = l[1:]
		if key.find(lookfor) > -1:
			key = trimto
		continue
	if not trimto:
		trimto = key

	if key in sequences:
		sequences[key] = sequences[key] + l
	else:
		sequences[key] = l


keep = []
	
for i in range(len(sequences[trimto])):
	if sequences[trimto][i] == '-':
		continue
	else:
		keep.append(i)

sequence = ""
key = 'n/a'
for l in infile:
	if l.find('>') == 0:
		seq = ""
		if key == 'n/a':
			key = l[1:]
			sequence = ""
			continue
		for i in range(len(keep)):
			seq = seq + sequence[keep[i]]
		print seq
		key = l[1:]
		sequence = ""
		continue
	sequence = sequence + l.replace('X', '-')

seq = ''
for i in range(len(keep)):
	seq = seq + sequence[keep[i]]
print seq
