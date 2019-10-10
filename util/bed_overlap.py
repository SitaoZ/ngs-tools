import argparse

def bed_overlap(bed):
    Bed = {}
    Count = {}
    with open(bed,'r') as F:
        for line in F:
            tmp = line.strip().split()
            chrid,start,end = tmp[0],tmp[1],tmp[2]
            key = "|".join([chrid,start,end])
            Bed[key]=[int(start),int(end)]
    for k1 in Bed.keys():
        s1,e1 = Bed[k1]
        for k2 in Bed.keys():
            s2,e2=Bed[k2]
            if e1 < s2 or e2 < s1:
                pass
            else:
                if Count.get(k1,0):
                    Count[k1] += 1
                else:
                    Count[k1] = 1
    for k in sorted(Count,key = lambda k:Count[k],reverse=True):
        out = k.split('|')
        out.append(str(Count[k]-1))
        print ('\t'.join(out))
        break
        
            
def main():
    """
    %prog [-options]
    The program for seq count analysis
    :return:
    """
    parser = argparse.ArgumentParser(prog='seq_count.py')
    parser.add_argument('-i', '--inputbed', help='the in file')
    args = parser.parse_args()
    seq_in = args.inputbed
    bed_overlap(seq_in)

if __name__=="__main__":
    main()
