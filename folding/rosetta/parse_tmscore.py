def read(scorefile):

    result = {}

    for line in scorefile:
        if line.startswith('RMSD'):
            result['RMSD'] = float(line.strip().split()[-1])
        if line.startswith('TM-score'):
            result['TM-score'] = float(line.strip().split()[2])
        if line.startswith('MaxSub'):
            result['MaxSub'] = float(line.strip().split()[1])
        if line.startswith('GDT-TS'):
            result['GDT-TS'] = float(line.strip().split()[1])
        if line.startswith('GDT-HA'):
            result['GDT-HA'] = float(line.strip().split()[1])

    scorefile.close()

    return result
