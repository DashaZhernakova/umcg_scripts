import subprocess
import sys

with open(sys.argv[1]) as f:
    header = f.readline()
    print(header + "\tzscores_per_cohort\tpvalues_per_cohort\theterogeneity_pvalue")
    colnames = {col : i for i,col in enumerate(header.rstrip().split("\t"))}
    snp_col = colnames["SNPCoord"]
    sv_col = colnames["ProbeName"]
    for l in f:
        spl = l.rstrip().split("\t")
        snp = spl[snp_col]
        sv = spl[sv_col]
        output = subprocess.check_output(['sh', '~/scripts/umcg_scripts/SV_GWAS/lookup_dSV.sh', sv, snp])
        res = output.split("\n")[-2].split(" ")
        spl.extend(res)
        print("\t".join(spl))
