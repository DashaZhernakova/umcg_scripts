#!/usr/bin/local/python
import sys
import gzip

in_fname = sys.argv[1]
rename_fname = sys.argv[2]
col_old = 0
col_new = 1
if len(sys.argv) > 3:
    col_old = int(sys.argv[3])
    col_new = int(sys.argv[4])

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

out_str = ""
header = in_file.readline().rstrip().split("\t")
for h in header:
    h2 = h
    if len(h) > 0:
        if h in rename_dict.keys():
             h2 = rename_dict[h]
        else:
             sys.stderr.write("WARNING! No such value in the renaming file: " + h + ". The old column name will be used!\n")
    out_str += "\t" + h2
print(out_str.replace("\t", "", 1))

for l in in_file:
    print(l.rstrip())
