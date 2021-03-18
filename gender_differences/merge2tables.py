import sys
from collections import defaultdict

fname1 = sys.argv[1]
fname2 = sys.argv[2]

n_res_cols = int(sys.argv[3])

all_samples = set()
dict_all = {}
with open(fname1) as f:
    header1 = f.readline().rstrip().split("\t")
    ncols1 = len(header1) - 1
    ncols2 = n_res_cols - ncols1
    for l in f:
        spl = l.rstrip().split("\t")
        sid = spl[0]
        all_samples.add(sid)
        dict_all[sid] = spl[1:] + ['NA']*ncols2

with open(fname2) as f:
    header2 = f.readline().rstrip().split("\t")
    assert(len(header2) - 1 == ncols2)
    for l in f:
        spl = l.rstrip().split("\t")
        sid = spl[0]
        if sid in all_samples:
            cur_entry = dict_all[sid]
            cur_entry[ncols1:] = spl[1:]
            dict_all[sid] = cur_entry
        else:
            dict_all[sid] = ['NA']*ncols1 + spl[1:]

for sid, dict_entry in dict_all.items():
    print (sid + "\t" + "\t".join(dict_entry))