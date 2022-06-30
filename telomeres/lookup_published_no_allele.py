import gzip
import sys
from collections import defaultdict

res = defaultdict(dict)
mtl_order = ["MTL_lymp", "MTL_gran", "MTL_CD45+_CD20-", "MTL_CD57+", "MTL_CD20+", "MTL_CD45-"]

publ = set()
with open(sys.argv[1], 'rt') as publ_f:
    for l in publ_f:
        spl = l.rstrip().split("\t")
        publ.add(l.strip())

allele_dict = {}
fname = "eQTLs.txt.gz"
with gzip.open(fname, 'rt') as f:
    f.readline()
    for l in f:
        spl = l.rstrip().split("\t")
        snp = spl[1]
        if snp in publ:
            mtl = spl[4]
            alleles = set(spl[8].split("/"))
            ea = spl[9]
            oa = alleles.difference(ea).pop()
            allele_dict[snp] = ea + "/" + oa
            beta = spl[17]
            res[snp][mtl] = (spl[0], beta)
print("\t".join(mtl_order))
for snp, qtls in res.items():
    out_line = snp + "\t" + allele_dict[snp]
    for mtl in mtl_order:
        cur_res = qtls[mtl]
        out_line += "\t" + cur_res[0] + "\t" + cur_res[1]
    print(out_line)
