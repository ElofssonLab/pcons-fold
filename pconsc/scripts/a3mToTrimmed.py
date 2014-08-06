#!/usr/bin/env python

import sys, os

infilef = sys.argv[1]
infile = open(infilef)

counter = 0
for l in infile:
    if '>' in l and not counter == 0:
        sys.stdout.write('\n>sequence{0:07d}/1-100\n'.format(counter))
        #sys.stdout.write('>sequence{0:07d}/1-100\n'.format(counter))
        counter += 1
    elif '>' not in l:
        l = l.strip()
        upperseq = ''.join([c for c in l if not c.islower()])
        upperseq = upperseq.replace('X', '-')
        sys.stdout.write(upperseq)
    elif '>' in l and counter == 0:
        sys.stdout.write('>target/1-100\n')
        counter += 1

sys.stdout.write('\n')
