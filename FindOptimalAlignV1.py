# This is for finding optimal alignment between 2 sequences.

import numpy as np

# Scoring Scheme
match = 2
mismatch = 0
gap = -2

def mm(a, b) :
    """Check a and b are match or mismatch"""
    if a == b :
        return match
    else :
        return mismatch

def find_optical_align_v1(sq1, sq2) :
    """print matching score and optimal alignment"""
    if len(sq2) > len(sq1) :
        sq1, sq2 = sq2, sq1
    sq1 = list(sq1)
    sq2 = list(sq2)
    sq1.insert(0, 0)
    sq2.insert(0, 0)

    mtx = np.zeros((len(sq1), len(sq2)))

    for i in range(0, len(sq1)) :
        mtx[i][0] = i * gap


    for i in range(0, len(sq2)) :
        mtx[0][i] = i * gap

    r, c = 1, 1
    while r < len(sq1) :
        while c < len(sq2) :
            mtx[r][c] = max([mtx[r - 1][c - 1] + mm(sq1[r], sq2[c]), mtx[r - 1][c] + gap, mtx[r][c - 1] + gap])
            c += 1
        r += 1
        c = 1

    r, c = len(sq1) - 1, len(sq2) - 1
    ssq1 = []
    ssq2 = []
    while (r, c) != (0, 0) :
        if mtx[r][c] == mtx[r - 1][c - 1] + mm(sq1[r], sq2[c]) :
            ssq1.append(sq1[r])
            ssq2.append(sq2[c])
            r, c = r - 1, c - 1
        elif mtx[r][c] == mtx[r - 1][c] + gap :
            ssq1.append(sq1[r])
            ssq2.append('-')
            r = r - 1
        else :
            ssq1.append('-')
            ssq2.append(sq2[c])
            c = c - 1

    ssq1 = ssq1[::-1]
    ssq2 = ssq2[::-1]

    print("*" * 50)
    print(f'[1] + Matching Score : {int(mtx[len(sq1) - 1][len(sq2) - 1])}')
    print('\n')
    print("*" * 50)
    print(f'[2] + One of the optimal Alignments')
    print(f'    {ssq1}')
    print(f'    {ssq2}')
    print('\n')

sq1 = "AAGCGGTGA"
sq2 = "AAGCGG"
find_optical_align_v1(sq2, sq1)
