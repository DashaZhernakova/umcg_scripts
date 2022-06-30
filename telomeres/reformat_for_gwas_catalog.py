import sys
import gzip
import re


out_path = sys.argv[2]
phenos = {"MTL_lymp" : "lymphocytes", "MTL_gran" : "granulocytes",  "MTL_CD45+_CD20-" : "naive_T-cells", "MTL_CD45-" : "memory_T-cells", "MTL_CD20+" : "B-cells", "MTL_CD57+" : "NK-cells"}
out_files = {}
header = "variant_id\tp_value\tchromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\n"


for k,v in phenos.items():
    out_f = open(out_path + "/" + v + ".summary_statistics.tsv", "w")
    out_files[k] = out_f
    out_f.write(header)


cnt = 1
with gzip.open(sys.argv[1]) as f:
    f.readline()
    for l in f:
        spl = l.rstrip().split("\t")
        pheno = spl[4]
        rsid = spl[1]
        if not spl[1].startswith("rs"):
            rsid = "NA"
        alleles = set(spl[8].split("/"))
        ea = spl[9]
        oa = alleles - set(ea)
        beta_se = spl[19]
        m = re.search(r"([\-\.0-9]+) \(([\-\.0-9]+)\)", beta_se)
        res_line = "\t".join([rsid, spl[0], spl[2], spl[3], ea, oa.pop(), m.group(1), m.group(2)])
        out_files[pheno].write(res_line + "\n")

        #cnt += 1
        #if cnt > 20:
        #    break

for k,v in out_files.items():
    out_files[k].close()
