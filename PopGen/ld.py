def _calc_burrows_pair(indp, indq):
    """Calculate burrows between 2 loci

    >> _calc_burrows_pair([(0, 1)], [(0, 1)])
    """
    fp = {}
    fq = {}
    Pijij = {}
    PijiJ = {}
    PijIj = {}
    PijIJ = {}
    for a1, a2 in indp:
        fp[a1] = fp.get(a1, 0) + 1
        fp[a2] = fp.get(a2, 0) + 1
    fp = {k: v / (2 * len(indp)) for k, v in fp.items()}
    for a1, a2 in indq:
        fq[a1] = fq.get(a1, 0) + 1
        fq[a2] = fq.get(a2, 0) + 1
    fq = {k: v / (2 * len(indq)) for k, v in fq.items()}
    print (fp)
    print(fq)
    for p in fp:
        for q in fq:
            Pijij[p, q] = 0
            PijiJ[p, q] = 0
            PijIj[p, q] = 0
            PijIJ[p, q] = 0
    for i in range(len(indp)):
        p0 = indp[i][0]
        p1 = indp[i][1]
        q0 = indq[i][0]
        q1 = indq[i][1]
        if p0 == p1 and q0 == q1:
            Pijij[p0, q0] = Pijij.get((p0, q0), 0) + 1
        elif p0 == p1:
            PijiJ[p0, q0] = PijiJ.get((p0, q0), 0) + 1
            PijiJ[p0, q1] = PijiJ.get((p0, q1), 0) + 1
        elif q0 == q1:
            PijIj[p0, q0] = PijIj.get((p0, q0), 0) + 1
            PijIj[p1, q0] = PijIj.get((p1, q0), 0) + 1
        else:
            PijIJ[p0, q0] = PijIJ.get((p0, q0), 0) + 1
            PijIJ[p0, q1] = PijIJ.get((p0, q1), 0) + 1
            PijIJ[p1, q0] = PijIJ.get((p1, q0), 0) + 1
            PijIJ[p1, q1] = PijIJ.get((p1, q1), 0) + 1
    Pijij = {k: 4 * v / len(indp) for k, v in Pijij.items()}
    PijiJ = {k: 2 * v / len(indp) for k, v in PijiJ.items()}
    PijIj = {k: 2 * v / len(indp) for k, v in PijIj.items()}
    PijIJ = {k: v / len(indp) for k, v in PijIJ.items()}
    p = list(fp.keys())[0]
    q = list(fq.keys())[0]
    rs = 0
    for p in fp:
        for q in fq:
            rs += Pijij[p, q] + PijIj[p, q] + PijiJ[p, q] + PijIJ[p, q] -\
                2 * fp[p] * fq[q]
    return rs
if __name__=="__main__":
    print(_calc_burrows_pair([(0, 1)], [(0, 1)]))