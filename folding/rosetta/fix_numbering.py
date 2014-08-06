import sys

sys.path.append('/bubo/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')
from Bio import pairwise2

import parse_pdb



def fix(pdb1_filename, pdb2_filename, out_filename, chain1='', chain2=''):

    if not chain1:
        chain1 = parse_pdb.get_first_chain(open(pdb1_filename, 'r'))
    pdb1 = parse_pdb.read(open(pdb1_filename, 'r'), chain1)
    seq1 = parse_pdb.get_atom_seq(open(pdb1_filename, 'r'), chain1)

    if not chain2:
        chain2 = parse_pdb.get_first_chain(open(pdb2_filename, 'r'))
    pdb2 = parse_pdb.read(open(pdb2_filename, 'r'), chain2)
    seq2 = parse_pdb.get_atom_seq(open(pdb2_filename, 'r'), chain2)

    align = pairwise2.align.globalms(seq1, seq2, 2, -1, -0.5, -0.1)
    #print seq1
    #print seq2

    #print align
    seq1_ali = align[-1][0]
    seq2_ali = align[-1][1]

    #print pdb2
    pdb2_idx = []
    offset = 0
    for i in xrange(len(seq2_ali)):
        if seq1_ali[i] == '-':
            offset -= 1
            idx = i+1 + offset
            pdb2_idx.append(idx)
        elif seq2_ali[i] == '-':
            continue
            #offset += 1
            #idx = i+1 + offset
            #pdb2_idx.append(idx)
        else:
            idx = i+1 + offset
            pdb2_idx.append(idx)
        #else:
    pdb2_new = ['', [], pdb2[2]]
    i = 0
    prev_idx = -1
    #print len(pdb2_idx)
    #print len(pdb2[1])
    for res in pdb2[1]:
        if i >= len(pdb2_idx):
            break
        new_res = []
        new_idx = pdb2_idx[i]
        if new_idx == 0:
            i = i+1
            continue
        elif new_idx == prev_idx:
            i = i+1
            continue
        else:
            for atm in res:
                new_idx_str = str(pdb2_idx[i])
                #print atm
                #print new_idx_str
                lendiff = 4 - len(new_idx_str)
                new_atm = atm[:22] + lendiff * ' ' + new_idx_str + atm[26:]
                new_res.append(new_atm)
            pdb2_new[1].append(new_res)

        prev_idx = new_idx
        i = i+1

    #print pdb1_filename
    #print pdb2_filename
    #print pdb2_idx
    #print len(pdb2_idx)
    #print align[-1]
    #print len(align[-1][1])
    if out_filename:
        pdb2_outfile = open(out_filename, 'w')
    else:
        pdb2_outfile = open('.'.join(pdb2_filename.split('.')[:-1]) + '.aligned.pdb', 'w')
    #print pdb2_new
    parse_pdb.write(pdb2_new, pdb2_outfile)
    


if __name__ == '__main__':

    pdb1_filename = sys.argv[1]
    pdb2_filename = sys.argv[2]
    outfilename = sys.argv[3]

    fix(pdb1_filename, pdb2_filename, outfilename)
