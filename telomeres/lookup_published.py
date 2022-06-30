import gzip
import sys
from collections import defaultdict

mtl_order = ["MTL_lymp", "MTL_gran", "MTL_CD45+_CD20-", "MTL_CD57+", "MTL_CD20+", "MTL_CD45-"]

#publ = {"rs3219104" : "A", "rs228596" : "C", "rs7776744" : "G", "rs10936599" : "T", "rs2736100" : "A", "rs7675998" : "A", "rs9420907" : "A", "rs8105767" : "A", "rs755017" : "A", "rs11125529" : "C", "rs412658" : "T", "rs3027234" : "T", "rs2535913" : "A"}
#publ = {"rs3219104" : "C", "rs10936600" : "T", "rs4691895" : "C", "rs7705526" : "A", "rs2853677" : "A", "rs59294613" : "A", "rs9419958" : "C", "rs228595" : "A", "rs2302588" : "C", "rs7194734" : "T", "rs8105767" : "G", "rs75691080" : "T", "rs34978822" : "G", "rs73624724" : "C", "rs55749605" : "A", "rs13137667" : "C", "rs34991172" : "G", "rs2736176" : "C", "rs3785074" : "G", "rs62053580" : "G"}

publ = {}
with open(sys.argv[1], 'rt') as publ_f:
    for l in publ_f:
        spl = l.rstrip().split("\t")
        publ[spl[0]] = spl[1]

res = defaultdict(dict)
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
            if publ[snp] == ea:
                beta = spl[17]
            elif publ[snp] == oa:
                beta = str(-1*float(spl[17]))
            else:
                print ("alleles are different for", snp, ": ", alleles, publ[snp])
            res[snp][mtl] = (spl[0], beta)

print("\t".join(mtl_order))
for snp, qtls in res.items():
    out_line = snp
    for mtl in mtl_order:
        cur_res = qtls[mtl]
        out_line += "\t" + cur_res[0] + "\t" + cur_res[1]
    print(out_line)
