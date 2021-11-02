#!/usr/bin/local/python
import sys

in_fname = sys.argv[1]
rename_fname = sys.argv[2]
col_old = 0
col_new = 1
if len(sys.argv) > 3:
    col_old = sys.argv[3]
    col_new = sys.argv[4]

rename_dict = {}
with open (rename_fname) as f:
    for l in f:
        spl = l.rstrip().split("\t")
        rename_dict[spl[col_old]] = spl[col_new]

with open(in_fname) as f:
    header = f.readline().rstrip().split("\t")
    for h in header:
        if len(h) > 0:
            h2 = rename_dict[h]
        else:
            h2 = h
        out_str += "\t" + h2
    print(out_str.replace("\t", "", 1))

    for l in f:
        print(l.rstrip())
