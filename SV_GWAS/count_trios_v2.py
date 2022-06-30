from collections import defaultdict


ids = {}
with open("ugli_to_dag3_ids.txt") as f:
    for l in f:
        spl = l.rstrip().split("\t")
        ids[spl[0]] = spl[1]


child_dict2 = {}
with open("/groups/umcg-lifelines/tmp01/projects/ov20_0051/umcg-dzhernakova/LL_genotypes/UGLI_mgs_v1/UGLI_family_relationships.txt") as f:
    f.readline()
    for l in f:
        spl = l.rstrip().split("\t")
        ch = ids.get(spl[0], "NA")
        f = ids.get(spl[1], "NA")
        m = ids.get(spl[2], "NA")
        if not ch == "NA" and not (f == "NA" and m == "NA"):
            child_dict2[ch] = [m,f]

child_dict3 = {}
with open("/groups/umcg-lifelines/tmp01/projects/ov20_0051/umcg-dzhernakova/LL_genotypes/UGLI_mgs_v1/UGLI_family_relationships.txt") as f:
    h=f.readline()
    for l in f:
        spl = l.rstrip().split("\t")
        if not spl[0] == "NA" and not (spl[1] == "NA" and spl[2] == "NA"):
            child_dict3[spl[0]] = spl[1:]


ids_ll = {}
#with open("/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file_project_pseudo_id.txt") as f:
with open("/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat") as f:
    for l in f:
        spl = l.rstrip().split("\t")
        ids_ll[spl[0]] = spl[1]

ids_dag3 = {}
with open("DAG3_link_microbiome") as f:
    for l in f:
        spl = l.rstrip().split("\t")
        ids_dag3[spl[1]] = spl[0]

child_dict3 = defaultdict(set)
#with open("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/family_relations/sec_family_relations.txt") as f:
out_f = open("DAG3_families_280622.txt", "w")
with open("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-rgacesa/DAG3_data_ready/family_data/LifeLines_families_info_withspouse.dat") as f:
    h = f.readline()
    for l in f:
        spl = l.rstrip().split("\t")
        #ch = ids.get(ids_ll.get(spl[1], "NA"), "NA")
        if (spl[2] != "0" and spl[3] != "0"):
            ch = ids_dag3.get(spl[0], "NA")
            f = ids_dag3.get(spl[2], "NA")
            m = ids_dag3.get(spl[3], "NA")
            out_f.write(ch +"\t" + spl[1] + "\t" +  f + "\t" + m + "\n")
            
            if not ch == "NA":
                child_dict3[spl[1]].add(ch)
                if not f == "NA":
                    child_dict3[spl[1]].add(f)
                if not m == "NA":
                    child_dict3[spl[1]].add(m)


out_f.close()

cnt = 1
for k,v in child_dict3.items():
    if len(v) > 2:
        cnt += 1
        print(k,v)
    if cnt > 5:
        break
cnt




