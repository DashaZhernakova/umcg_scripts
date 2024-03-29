import sys
from collections import defaultdict

fname1 = sys.argv[1]
fname2 = sys.argv[2]

with open(fname1) as f:
    header1 = f.readline().rstrip().split("\t")
    ncols1 = len(header1) - 1
with open(fname2) as f:
    header2 = f.readline().rstrip().split("\t")
    ncols2 = len(header2) - 1
sys.stderr.write("Numberr of columns in " + fname1 + ": " + str(ncols1) + "\n")
sys.stderr.write("Numberr of columns in " + fname2 + ": " + str(ncols2) + "\n")

all_samples = set()
dict_all = {}
with open(fname1) as f:
    f.readline()
    for l in f:
        spl = l.rstrip().split("\t")
        sid = spl[0]
        all_samples.add(sid)
        dict_all[sid] = spl[1:] + ['NA']*ncols2

with open(fname2) as f:
    f.readline()
    for l in f:
        spl = l.rstrip().split("\t")
        sid = spl[0]
        if sid in all_samples:
            cur_entry = dict_all[sid]
            cur_entry[ncols1:] = spl[1:]
            dict_all[sid] = cur_entry
        else:
            dict_all[sid] = ['NA']*ncols1 + spl[1:]

print("\t".join(header1) + "\t" + "\t".join(header2[1:]))
for sid, dict_entry in dict_all.items():
    print (sid + "\t" + "\t".join(dict_entry))