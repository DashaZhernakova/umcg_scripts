#!/usr/local/bin/python
from collections import defaultdict

def comparator(x, y):
    #ss_col = 3
    #nsnp_col = 5
    ss_col = 4
    nsnp_col = 6
    if int(x[ss_col]) > int(y[ss_col]):
        return -1
    if int(x[ss_col]) < int(y[ss_col]):
        return 1
    if int(x[ss_col]) == int(y[ss_col]):
        if int(x[nsnp_col]) > int(y[nsnp_col]):
            return -1
        if int(x[nsnp_col]) < int(y[nsnp_col]):
            return 1
        if int(x[nsnp_col]) == int(y[nsnp_col]):
            return 0

out = open("/Users/dashazhernakova/Documents/UMCG/data/MR/results2/mibiogen/mibiogenSep2019/tmp_res4filtering_res.txt", "w")
pair2line = defaultdict(list)
with open("/Users/dashazhernakova/Documents/UMCG/data/MR/results2/mibiogen/mibiogenSep2019/tmp_res4filtering.txt") as f:
    out.write(f.readline().strip() + "\tfilter_samplesize\n")
    for l in f:
        spl = l.strip().split("\t")
        spl.append("FALSE")
        pair2line[spl[0] + ":" + spl[1]].append(spl)


#print pair2line
for p, spl_list in pair2line.items():
    spl_list.sort(comparator)
    spl_list[0][-1] = "TRUE"
    for sp in spl_list:
        out.write("\t".join(sp) + "\n")
out.close()