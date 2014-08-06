#!/usr/bin/env python

import sys, os

infilef = sys.argv[1]
infile = open(infilef)

count = 0
line = ""

for l in infile:
    if '>' in l:
        if (count != 0):
            sys.stdout.write(line)
            sys.stdout.write('\n')
            line = ""

        count += 1
        continue

    l = l.strip()
    
    upperseq = ''.join([c for c in l if not c.islower()])
    upperseq = upperseq.replace('X', '-')
    line = line + upperseq

sys.stdout.write(line)
sys.stdout.write('\n')
