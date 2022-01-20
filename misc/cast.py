#!/usr/bin/python
import sys
from collections import defaultdict


if sys.argv[1] != "-":
    f = open(sys.argv[1])
else:
    f = sys.stdin

col_key = sys.argv[2]
col_val = sys.argv[3]

map_dict = defaultdict(set)
for line in f:
	spl = line.rstrip().split("\t")
	map_dict[spl[col_key]].add(spl[col_val])


for k, v in map_dict.items():
    print (k + "\t" + ",".join(v))