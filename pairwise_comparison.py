import pandas as pd
import numpy as np
import os

def comp(smp1, smp2) :
    count = 0
    n1 = len(smp1)
    lst1 = []
    for i in positions1 :
        lst1.append(smp1['gene'][i] + str(smp1['chr'][i])+ str(smp1['position'][i])+ str(smp1['ref'][i])+ str(smp1['alt'][i]) + str(smp1['variant'][i]))

    n2 = len(smp2)
    lst2 = []
    for i in positions2 :
        lst2.append(smp2['gene'][i] + str(smp2['chr'][i])+ str(smp2['position'][i])+ str(smp2['ref'][i])+ str(smp2['alt'][i]) + str(smp2['variant'][i]))

    for i in range(0, n1) :
        if lst1[i] in lst2 :
            count += 1
    return count

adr = "./nsyn_2E-05.xlsx"
data = pd.read_excel(adr)
data_np = pd.DataFrame.to_numpy(data)
data_np = pd.DataFrame(data_np, columns = ["gnomAD", "gene", "smp", "chr", "position", "ref", "alt", "variant"])
sample = []
for i in range(len(data_np)) :
    if data_np["smp"][i] not in sample :
        sample.append(data_np["smp"][i])
sample.sort()
sample.insert(0, "N")
mtx = np.zeros((len(sample),len(sample)))
mtx = pd.DataFrame(mtx, index=sample, columns=sample)
sample.pop(0)

dup = []
for a in sample :
    for b in sample :
        if a == b or b in dup :
            mtx[b][a] = "."
        else :
            data1 = data_np['smp'] == a
            positions1 = np.flatnonzero(data1)
            tmp1 = data_np.iloc[positions1]
            data2 = data_np['smp'] == b
            positions2 = np.flatnonzero(data2)
            tmp2 = data_np.iloc[positions2]
            rst = comp(tmp1, tmp2)
            mtx[a]["N"] = len(tmp1)
            mtx["N"][a] = len(tmp1)
            mtx[b]["N"] = len(tmp2)
            mtx["N"][b] = len(tmp2)
            tmp = (rst/(len(tmp1) + 1) + rst/(len(tmp2) + 1))/2
            mtx[b][a] = round(tmp, 2)
            dup.append(a)

mtx.to_excel("result.xlsx")