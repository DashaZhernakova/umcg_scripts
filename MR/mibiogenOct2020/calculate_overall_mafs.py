import sys
import os
import gzip

def get_AF(spl, mr_alleles):
    maf_alleles = spl[2].split("/")
    AC = 0
    AN = 0
    geno_cnts = [ int(spl[3].split(" ")[0]), int(spl[4].split(" ")[0]), int(spl[5].split(" ")[0]) ]
    if maf_alleles == mr_alleles:
        AN = 2*sum(geno_cnts)
        AC = geno_cnts[1] + 2*geno_cnts[2]
    elif (maf_alleles[1] == mr_alleles[0] and maf_alleles[0] == mr_alleles[1]):
        AN = 2*sum(geno_cnts)
        AC = geno_cnts[1] + 2*geno_cnts[0]
    else:
        print("wrong alleles!", spl[0], mr_alleles, maf_alleles)
    return(AC, AN)


mr_dir = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/data/mibiogenOct2020/"
maf_dir = "/groups/umcg-lifelines/tmp03/users/umcg-akurilshchikov/MiBioGen/MAFs/"
res_fname = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/data/mibiogenOct2020/all_snp_afs.txt"
#mr_dir = "C:/Users/Dasha/work/UMCG/data/MR/results2/mibiogen/mibiogenOct2020/test/res/"
#maf_dir = "C:/Users/Dasha/work/UMCG/data/MR/results2/mibiogen/mibiogenOct2020/test/mafs/"
mr_snp_dict = {}

for mr_fname in os.listdir(mr_dir):
    with gzip.open(mr_dir + mr_fname, 'rt') as mr_file:
        mr_file.readline()
        for l in mr_file:
            spl = l.rstrip().split("\t")
            mr_snp_dict.setdefault(spl[1] + ":" + spl[2], [spl[4] + "/" + spl[5], 0, 0])
print("Read", len(mr_snp_dict), "SNPs")

for maf_fname in os.listdir(maf_dir):
    with gzip.open(maf_dir + maf_fname, 'rt') as maf_file:
        maf_file.readline()
        maf_file.readline()
        for l in maf_file:
            spl = l.rstrip().split("\t")
            mr_snp = mr_snp_dict.get(spl[0])
            if not mr_snp or spl[9] == "false":
                continue
            mr_alleles = mr_snp[0].split("/")
            AC, AN = get_AF(spl, mr_alleles)
            mr_snp[1] += AC
            mr_snp[2] += AN
            mr_snp_dict[spl[0]] = mr_snp

res_file = open(res_fname, "w")
res_file.write("SNP\tother_allele\teffect_allele\teaf\n")
for snp, lst in mr_snp_dict.items():
    res_file.write(snp + "\t" + "\t".join(lst[0].split("/")) + "\t" + str(lst[1]/lst[2]) + "\n")
res_file.close()

for mr_fname in os.listdir(mr_dir):
    with gzip.open(mr_dir + mr_fname, 'rt') as mr_file:
        res_mr_file = open(mr_dir + mr_fname.replace("summary.txt.gz", "summary.afs.txt"), "w")
        l = mr_file.readline()
        res_mr_file.write(l.rstrip() + "\teaf\n")
        for l in mr_file:
            spl = l.rstrip().split("\t")
            mr_snp = mr_snp_dict.get(spl[1] + ":" + spl[2])
            af = str(mr_snp[1]/mr_snp[2])
            res_mr_file.write(l.rstrip() + "\t" + af + "\n")
        res_mr_file.close()
print("Finished")