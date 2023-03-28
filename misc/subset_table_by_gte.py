import sys
import gzip
import pandas as pd

table_fname = sys.argv[1]
gte_fname = sys.argv[2]
out_fname = sys.argv[3]
if len(sys.argv) > 4:
    col = int(sys.argv[4])
else:
    col = 1

table = pd.read_csv(table_fname, sep = "\t", index_col = 0)
gte_f = open(gte_fname, "rt")

print("Original number of rows and columns: {}.".format(table.shape))

samples = [l.rstrip().split("\t")[col] for l in gte_f.readlines()]
samples_overlap = [x for x in samples if x in table.columns]
table_subset = table.loc[:,samples_overlap]

print("New number of rows and columns: {}.".format(table_subset.shape))

table_subset.to_csv(out_fname, sep="\t")

