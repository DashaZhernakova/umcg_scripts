#!/usr/bin/python
import sys
from collections import defaultdict


if sys.argv[1] != "-":
    f = open(sys.argv[1])
else:
    f = sys.stdin

col_key = int(sys.argv[2])
col_val = int(sys.argv[3])

col_val2 = None
if (len(sys.argv) > 4):
    col_val2 = int(sys.argv[4])
    dict2 = defaultdict(set)

map_dict = defaultdict(set)
for line in f:
    spl = line.rstrip().split("\t")
    map_dict[spl[col_key]].add(spl[col_val])
    if col_val2 is not None:
        dict2[spl[col_key]].add(spl[col_val2])


for k, v in map_dict.items():
    if col_val2 is None:
        print (k + "\t" + ",".join(v))
    else:
        second = dict2.get(k, "NA")
        print (k + "\t" + ",".join(v) + "\t" + ",".join(second))