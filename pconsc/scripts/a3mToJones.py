#!/usr/bin/env python

import sys, os

infilef = sys.argv[1]
infile = open(infilef)

for l in infile:
    if '>' in l:
        continue
    upperseq = ''.join([c for c in l if not c.islower()])
    upperseq = upperseq.replace('X', '-')
    sys.stdout.write(upperseq)
