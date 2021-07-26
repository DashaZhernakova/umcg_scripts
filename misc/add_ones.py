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


if fname1 != "stdin":
    if fname1.endswith(".gz"):
        i_file = gzip.open(fname1)
    else:
        i_file = open(fname1)
else:
    i_file = sys.stdin


if new_col_name:
    print (f.readline().rstrip() + "\t" + new_col_name)
for l in i_file:
    spl = l.rstrip().split("\t")
    add_val = '0'
    if spl[0] in id_set:
        add_val = '1'
    print (l.rstrip() + "\t" + add_val)
