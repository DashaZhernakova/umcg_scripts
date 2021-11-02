#!/usr/bin/local/python
import sys
import gzip

in_fname = sys.argv[1]
rename_fname = sys.argv[2]
col_old = 0
col_new = 1
if len(sys.argv) > 3:
    col_old = sys.argv[3]
    col_new = sys.argv[4]

if in_fname != "stdin":
    if in_fname.endswith(".gz"):
        in_file = gzip.open(in_fname)
    else:
        in_file = open(in_fname)
else:
    in_file = sys.stdin


rename_dict = {}
with open (rename_fname) as f:
    for l in f:
        spl = l.rstrip().split("\t")
        rename_dict[spl[col_old]] = spl[col_new]


header = in_file.readline().rstrip().split("\t")
for h in header:
    if len(h) > 0:
        h2 = rename_dict[h]
    else:
        h2 = h
    out_str += "\t" + h2
print(out_str.replace("\t", "", 1))

for l in in_file:
    print(l.rstrip())
