#!/usr/bin/env python

import sys, os, random

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

	if key in sequences:
		sequences[key] = sequences[key] + l
	else:
		sequences[key] = l


chosen = ['target']
for l in infile:
	if l.find('>') == 0:
	        key = l[1:]
                chosen.append(key)

keep = []
	
for i in range(len(sequences[trimto])):
	if sequences[trimto][i] == '-':
		continue
	else:
		keep.append(i)

maxseq = 1000000/len(keep)
if len(chosen) > maxseq:
        chosen2 = random.sample(chosen, maxseq)
        chosen2.append('target')
        chosen = chosen2

sequence = ""
key = 'n/a'
counter = 1
skip = False

for l in infile:
	if l.find('>') == 0:
		if key == 'n/a':
			seq = ""
			key = l[1:]
                	if key.find(lookfor) > -1:
                        	key = trimto
			sequence = ""
			continue
                
	
                if key in chosen:
			seq = ""
			if key.find('dssp') > -1:
				continue
		        for i in range(len(keep)):
		        	seq = seq + sequence[keep[i]]
                        if key != 'target':
                                key = 'sequence' + str(counter) + '/1-100'
                                counter = counter + 1
                        else:
                                key = 'target/1-100'
		        print '>' + key
        		print seq
	        key = l[1:]
		sequence = ""
                if key.find(lookfor) > -1:
                       	key = trimto
		continue

	sequence = sequence + l.replace('X', '-')

if key in chosen:
	seq = ""
	if key.find('dssp') < 0:
		for i in range(len(keep)):
			seq = seq + sequence[keep[i]]
		if key != 'target':
			key = 'sequence' + str(counter) + '/1-100'
			counter = counter + 1
		else:
			key = 'target/1-100'
		print '>' + key
		print seq
