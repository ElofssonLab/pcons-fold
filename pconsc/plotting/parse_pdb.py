import sys
import operator
import numpy as np
from collections import defaultdict


def parse_atm_record(line):

    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])
    
    return record


def read(pdbfile):

    header = ''
    res_lst = []
    atm_lst = []
    tail = ''

    seen_atoms = False
    curr_resi = 0
    prev_resi = 0
    
    for line in pdbfile:
        if not line.startswith('ATOM') and not seen_atoms:
            header += line
        elif not line.startswith('ATOM') and seen_atoms:
            tail += line
        else:
            atm_record = parse_atm_record(line)
            if not seen_atoms:
                curr_resi = atm_record['res_no']
                prev_resi = curr_resi
            seen_atoms = True
            curr_resi = atm_record['res_no']
            if curr_resi == prev_resi:
                atm_lst.append(line)
            else:
                #atm_lst.append(line)
                res_lst.append(atm_lst)
                atm_lst = [line]
            prev_resi = curr_resi
    res_lst.append(atm_lst)
     
    pdbfile.close()
    pdb_lst = [header, res_lst, tail]
    return pdb_lst


def write(pdb_lst, outfile):

    outfile.write(pdb_lst[0])

    for res in pdb_lst[1]:
        for atm in res:
            outfile.write(atm)
            
    outfile.write(pdb_lst[2])
    outfile.close()


def get_coordinates(pdbfile, chain):

    res_dict = defaultdict(list)

    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        if atm_record['chain'] != ' ' and atm_record['chain'] != chain:
            continue

        res_i = atm_record['res_no']
        atm = [atm_record['x'], atm_record['y'], atm_record['z']]

        """
        line_arr = line.split()
        #print line_arr

        if line_arr[2].startswith('H'):
            continue

        if len(line_arr[2]) > 4:
            if line_arr[3] != chain:
                continue
            try:
                res_i = int(line_arr[4])
            except ValueError as exc:
                continue
            try:
                atm = map(float, line_arr[5:8])
            except ValueError as exc:
                atm = [float('inf'), float('inf'), float('inf')]
        else:
            if line_arr[4] != chain:
                continue
            try:
                res_i = int(line_arr[5])
            except ValueError as exc:
                continue
            try:
                atm = map(float, line_arr[6:9])
            except ValueError as exc:
                atm = [float('inf'), float('inf'), float('inf')]
        """

        res_dict[res_i].append(np.array(atm))
        
    pdbfile.close()
    return sorted(res_dict.iteritems(), key=operator.itemgetter(0))


def get_cb_coordinates(pdbfile, chain):

    cb_lst = []
    res_dict = defaultdict(list)

    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)

        if atm_record['chain'] != ' ' and atm_record['chain'] != chain:
            continue

        res_i = atm_record['res_no']
        
        atm = [float('inf'), float('inf'), float('inf')]
        if 'GLY' in atm_record['res_name'] and atm_record['atm_name'] == 'CA':
                atm = [atm_record['x'], atm_record['y'], atm_record['z']]
                res_dict[res_i].append(np.array(atm))    
        elif atm_record['atm_name'] == 'CB':
                atm = [atm_record['x'], atm_record['y'], atm_record['z']]
                res_dict[res_i].append(np.array(atm))    
        

        """
        line_arr = line.split()

        if line_arr[4] != chain:
            continue

        if line_arr[2] == 'CA':
            try:
                res_i = int(line_arr[5])
            except ValueError as exc:
                continue
            try:
                atm = map(float, line_arr[6:9])
            except ValueError as exc:
                atm = [float('inf'), float('inf'), float('inf')]
            res_dict[res_i].append(np.array(atm))
        else:
            if line_arr[2] != 'CB':
                continue
            try:
                res_i = int(line_arr[5])
            except ValueError as exc:
                continue
            if len(line_arr[3]) > 3 and line_arr[3].startswith('A'):
                try:
                    atm = map(float, line_arr[6:9])
                except ValueError as exc:
                    atm = [float('inf'), float('inf'), float('inf')]
                atm_count += 1
                cb_lst.append(np.array(atm))
                res_dict[res_i].append(np.array(atm))
            elif len(line_arr[3]) == 3:
                try:
                    atm = map(float, line_arr[6:9])
                except ValueError as exc:
                    atm = [float('inf'), float('inf'), float('inf')]
                atm_count += 1
                cb_lst.append(np.array(atm))
                res_dict[res_i].append(np.array(atm))
        #print line_arr[2]
        
        if line_arr[3] == 'GLY' or line_arr[3] == 'AGLY':
            if line_arr[2] != 'CA':
                continue
            try:
                res_i = int(line_arr[5])
            except ValueError as exc:
                continue
            try:
                atm = map(float, line_arr[6:9])
            except ValueError as exc:
                atm = [float('inf'), float('inf'), float('inf')]
            #atm = map(float, line_arr[6:9])
            atm_count += 1
            cb_lst.append(np.array(atm))
        """
 
    cb_lst = []
    for i in xrange(res_i):
        if len(res_dict[i+1]) > 1:
            #print res_dict[i+1]
            cb_lst.append(res_dict[i+1][1])
        elif len(res_dict[i+1]) == 1:
            #print res_dict[i+1]
            cb_lst.append(res_dict[i+1][0])
    #print atm_count 
    pdbfile.close()
    #print cb_lst
    return cb_lst


def get_atom_seq(pdbfile, chain):

    three_to_one = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C', 'GLY':'G', 'PRO':'P', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 'UNK': 'X'}
    res_dict = {}
 
    if not chain:
        chain = get_first_chain(pdbfile)
        pdbfile.seek(0)

    res_name = ''
    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        if atm_record['chain'] != ' ' and atm_record['chain'] != chain:
            continue

        res_i = atm_record['res_no']
        #print atm_record['res_name']
        if atm_record['res_name'] in three_to_one:
            #res_name = three_to_one[atm_record['res_name']]
            #print res_name
            res_name = three_to_one[atm_record['res_name']]
        #else:
            #res_name = ''
            #continue

        res_dict[res_i] = res_name

    res_lst = sorted(res_dict.iteritems(), key=operator.itemgetter(0))
    atom_seq = ''

    for res in res_lst:
        atom_seq += res[1]

    pdbfile.close()
    return atom_seq


def get_first_chain(pdbfile):

    for line in pdbfile:
        if not line.startswith('ATOM'):
            continue
        atm_record = parse_atm_record(line)
        break

    return atm_record['chain']
 



if __name__ == '__main__':

    pdbfile = open(sys.argv[1], 'r')
    chain = sys.argv[2]
    #print get_atom_seq(pdbfile, chain)
    pdbfile.close()
    pdbfile = open(sys.argv[1], 'r')
    #print get_coordinates(pdbfile)
    pdbfile.close()
