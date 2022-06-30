
from collections import defaultdict



child_dict = {}
with open("/groups/umcg-lifelines/tmp01/projects/ov20_0051/umcg-dzhernakova/LL_genotypes/UGLI_mgs_v1/DAG3_family_relationships.txt") as f:
    f.readline()
    for l in f:
        spl = l.rstrip().split("\t")
        if not spl[0] == "NA" and not (spl[1] == "NA" and spl[2] == "NA"):
            child_dict[spl[0]] = spl[1:]



sv_dict = defaultdict(set)
with open("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/data/DAG3.dSV.filtered.txt") as f2:
        header = {i : sv for i, sv in enumerate(f2.readline().rstrip().split("\t"))}
        for l in f2:
            spl = l.rstrip().split("\t")
            for i, val in enumerate(spl):
                if not val == "NA":
                    sv_dict[i].add(spl[0])


for sv, sv_samples in sv_dict.items():
    sum_count = 0
    for child, parents in child_dict.items():
        count = 0
        if child in sv_samples:
            #sv_samples.remove(child)
            if parents[0] in sv_samples or parents[1] in sv_samples:
                count = 1
        sum_count += count
    print(header[sv] + "\t" + str(sum_count))


    
