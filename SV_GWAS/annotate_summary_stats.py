import gzip
from collections import defaultdict
import sys

svtype = "dSV"

annot_fname = "/data/umcg-tifn/SV/SV_GWAS/data" + svtype + "_annotation.txt"
annot_dict = {}
with open(annot_fname) as annot_f:
    annot_f.readline()
    for l in annot_f:
        spl = l.rstrip().split("\t")
        annot_dict[spl[8]] = [spl[2], spl[3]]

freq_dict = {}
cohorts = ["DAG3", "LLD", "300OB", "500FG"]
cohorts_dict = {name : num for num, name in enumerate(cohorts)}
for c_num, c in enumerate(cohorts):
    with gzip.open("/data/umcg-tifn/SV/SV_GWAS/genotypes/" + c + "/" + c + "_filtered.frq.counts.gz") as frq_f:
        frq_f.readline()
        cnt = 0
        for l in frq_f:
            spl = l.rstrip().split()
            spl[4] = float(spl[4])
            spl[5] = float(spl[5])
            af1 = spl[4]/(spl[4] + spl[5])
            af2 = spl[5]/(spl[4] + spl[5])
            res_dict = {spl[2] : af1, spl[3] : af2}
            if spl[1] in freq_dict:
                freq_dict[spl[1]][c_num] = res_dict
            else:
                res_lst = ["NA"]*4
                res_lst[c_num] = res_dict
                freq_dict[spl[1]] = res_lst

all_bacs_f = open("/data/umcg-tifn/SV/SV_GWAS/data/" + svtype + "_per_cohort.txt")
all_bacs_f.readline()
all_bacs = [l.split("\t")[0] for l in all_bacs_f]

for sv in all_bacs:
    #sv_annot = annot_dict[sv]

sv='B.longum:72'
with gzip.open("/data/umcg-tifn/SV/SV_GWAS/results_all_summary_stats/dSV/meta/" + sv + "/" + sv +".meta_res.annot.tbl.gz") as f:
    cnt = 0
    for l in f:
        spl = l.rstrip().split()
        #print(spl[1])
        freq_res = freq_dict[spl[1]]
        res_freq_lst = ["NA"]*4
        snp_cohorts = spl[8].split(",")
        cohorts_used = spl[7]
        for i, c in enumerate(snp_cohorts):
            if not cohorts_used[i] == '?':
                c_num = cohorts_dict[c]
                c_res = freq_res[c_num]
                if not c_res == 'NA':
                    res_freq_lst[c_num] = align_alleles(spl, c_res)
        print (l.rstrip() + "\t" + "\t".join(res_freq_lst))
        #cnt += 1
        #if cnt > 100:
        #    break

# think whether to output meta MAF
# replace DAG3 with DMP

def align_alleles(spl, c_res):
    ea = spl[2].upper()
    oa = spl[3].upper()
    if ea in c_res and oa in c_res:
        return(str(c_res[ea]))
    elif not isAT_GCsnp(ea,oa) and complement(ea) in c_res and complement(oa) in c_res:
        return(str(c_res[complement(ea)]))
    else:
        print ("\t".join(spl) + "\t" +"\t".join(c_res.keys()))
    return("NA")
    #elif ea not in c_res and oa not in c_res:
    #    if 


def complement(nucl):
    if nucl == 'A':
        return 'T'
    if nucl == 'T':
        return 'A'
    if nucl == 'G':
        return 'C'
    if nucl == 'C':
        return 'G'
 
def isAT_GCsnp(ref, alt):
    if ref == 'T' and alt == 'A':
        return True
    if ref == 'A' and alt == 'T':
        return True 
    if ref == 'G' and alt == 'C':
        return True
    if ref == 'C' and alt == 'G':
        return True
    return False


