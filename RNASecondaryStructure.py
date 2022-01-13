#  Finding maximum number of complementary base pairs in RNA sequence

import numpy as np

seqdic = {'A' : 'U', "U" : "A", "G" : "C", "C" : "G"}
def mm(a, b) :
    if seqdic[a] == b :
        return 1
    return 0

def nussinov_folding(seq1) :
    seq1 = list(seq1)
    mtx = np.zeros((len(seq1), len(seq1)))
    def bifurcation(r, c) :
        tmp = []
        for k in range(r + 1, c) :
            bi = mtx[r][k] + mtx[k + 1][c]
            tmp.append(bi)
        return max(tmp)

    i, j = 0, 1
    a = 2
    while True :
        if j == len(seq1) :
            i = 0
            j = a
            a += 1
        
        if j == len(seq1) :
            break

        if j - i < 2 :
            mtx[i][j] = max(mtx[i + 1][j], mtx[i][j - 1], mtx[i + 1][j - 1] + mm(seq1[i], seq1[j]))
        else :
            mtx[i][j] = max(mtx[i + 1][j], mtx[i][j - 1], mtx[i + 1][j - 1] + mm(seq1[i], seq1[j]), bifurcation(i, j))
        
        i += 1
        j += 1

    rcd = []
    stk = []
    stk.append((0, len(seq1) - 1))
    while stk :
        (i,j) = stk.pop()
        if i >= j :
            continue
        elif mtx[i+1][j] == mtx[i][j] :
            stk.append((i+1, j))
        elif mtx[i][j-1] == mtx[i][j] :
            stk.append((i, j-1))
        elif mtx[i+1][j-1] + mm(seq1[i], seq1[j]) :
            rcd.append((i, j))
            stk.append((i + 1, j - 1))
        else :
            for k in range(i + 1, j - 1) :
                if mtx[i][k] + mtx[k+1][j] == mtx[i][j] :
                    stk.append((k + 1, j))
                    stk.append((i, k))
                    break
    print(f"[+] Original Sequence : {''.join(seq1)}")
    print(f"[+] Maximum number of complementary base pairs : {int(mtx[0][len(seq1) - 1])}")
    print(f"[+] Expected linked base pairs : {rcd}")

nussinov_folding("GCGAGCG")