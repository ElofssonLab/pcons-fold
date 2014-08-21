import sys

def parse(f):
    ppv_dict = {}
    f.readline()
    for l in f:
        l_lst = l.strip().split()
        ppv_dict[l_lst[0]] = float(l_lst[1])
    return ppv_dict

if __name__ == "__main__":
    print parse(open(sys.argv[1],'r'))
