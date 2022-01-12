# This is for finding optimal alignment between 2 sequences.
# I added affine gap penalty from FindOptimalAlignV1

import numpy as np

# Scoring Scheme
match = 2
mismatch = 1
Ogap = -2
Egap = -1

def mm(a, b) :
    """Check a and b are match or mismatch"""
    if a == b :
        return match
    else :
        return mismatch

def find_optimal_align_v2(sq1, sq2) :
    """print matching score and optimal alignment"""
    if len(sq2) > len(sq1) :
        sq1, sq2 = sq2, sq1
    sq1 = list(sq1)
    sq2 = list(sq2)
    sq1.insert(0, 0)
    sq2.insert(0, 0)

    M = np.zeros((len(sq1), len(sq2)))
    IX = np.zeros((len(sq1), len(sq2)))
    IY = np.zeros((len(sq1), len(sq2)))

# -10000 means negative infinite
    for i in range(0, len(sq1)) :
        M[i][0] = -10000
        IY[i][0] = -10000
        if i == 0 :
            IX[i][0] = 0
        elif i == 1 :
            IX[i][0] = Ogap
        else :
            IX[i][0] = IX[i - 1][0] + Egap


    for i in range(0, len(sq2)) :
        M[0][i] = -10000
        IX[0][i] = -10000
        if i == 0 :
            IY[0][i] = 0
        elif i == 1 :
            IY[0][i] = Ogap
        else :
            IY[0][i] = IY[0][i - 1] + Egap

    M[0][0] = 0
    IX[0][0] = 0
    IY[0][0] = 0

    r, c = 1, 1
    while r < len(sq1) :
        while c < len(sq2) :
            M[r][c] = mm(sq1[r], sq2[c]) + max(M[r - 1][c - 1], IX[r - 1][c - 1], IY[r - 1][c - 1])
            IX[r][c] = max(M[r - 1][c] + Ogap, IX[r - 1][c] + Egap)
            IY[r][c] = max(M[r][c - 1] + Ogap, IY[r][c - 1] + Egap)
            c += 1
        r += 1
        c = 1

    r, c = len(sq1) - 1, len(sq2) - 1
    m = 0
    ix = 0
    iy = 0
    if M[r][c] == max(M[r][c], IX[r][c], IY[r][c]) :
        m  = 1
        ix = 0
        iy = 0
    elif IX[r][c] == max(M[r][c], IX[r][c], IY[r][c]) :
        m  = 0
        ix = 1
        iy = 0
    else :
        m  = 0
        ix = 0
        iy = 1

    ssq1 = []
    ssq2 = []
    while (r, c) != (0, 0) :
        if m == 1 :
            ssq1.append(sq1[r])
            ssq2.append(sq2[c])
            if M[r][c] == mm(sq1[r], sq2[c]) + IX[r - 1][c - 1] :
                m = 0
                ix = 1
            elif M[r][c] == mm(sq1[r], sq2[c]) + IY[r - 1][c - 1] :
                m = 0
                iy = 1
            r -= 1
            c -= 1

        elif ix == 1 :
            ssq1.append(sq1[r])
            ssq2.append('-')
            if IX[r][c] == M[r - 1][c] + Ogap :
                m = 1
                ix = 0
            r -= 1

        else :
            ssq1.append('-')
            ssq2.append(sq2[c])
            if IY[r][c] == M[r][c - 1] + Ogap :
                m = 1
                iy = 0
            c -=1
    

    ssq1 = ssq1[::-1]
    ssq2 = ssq2[::-1]

    print("*" * 50)
    print(f'[1] + Matching Score : {int(max(M[len(sq1) - 1][len(sq2) - 1], IX[len(sq1) - 1][len(sq2) - 1], IY[len(sq1) - 1][len(sq2) - 1]))}')
    print('\n')
    print("*" * 50)
    print(f'[2] + One of the optimal Alignments')
    print(f'    {ssq1}')
    print(f'    {ssq2}')
    print('\n')

sq1 = "AGT"
sq2 = "AAGC"
find_optimal_align_v2(sq2, sq1)