#!/usr/bin/python
import sys

fname1 = sys.argv[1]
fname2 = sys.argv[2]
if len(sys.argv) > 3:
    new_col_name = sys.argv[3]
else:
    new_col_name = None
id_set = set()
with open(fname2) as f:
    id_set = set(l.rstrip() for l in f)

with open(fname2) as f:
    if new_col_name:
        print (f.readline() + "\t" + new_col_name)
    for l in f:
        spl = l.rstrip().split("\t")
        add_val = 0
        if spl[0] in id_set:
            add_val = 1
        print (l.rstrip() + "\t" + add_val)