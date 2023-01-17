ml PLINK/1.9-beta6-20190617
ml KING

d=/groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/qc/
cd $d
script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#
# Mothers
#

# Relatedness IBD analysis

king -b ../all_chr.mothers.bed --related --prefix all_chr.mothers

grep -V "UN" all_chr.mothers.kin 
# FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    HetConc HomIBS0 Kinship IBD1Seg IBD2Seg PropIBD InfType Error
# 0       LLNEXT000145    LLNEXT005861    25096550        1.000   0.0000  0.0453  0.0030  0.4252  0.1265  0.2588  0.5218  0.2426  0.5035  FS      1
# 0       LLNEXT000307    LLNEXT134478    25096550        1.000   0.0000  0.0445  0.0033  0.4235  0.1366  0.2529  0.5100  0.2465  0.5015  FS      1
# 0       LLNEXT000316    LLNEXT201947    25096550        1.000   0.0000  0.0431  0.0027  0.3958  0.1139  0.2477  0.5950  0.1729  0.4704  FS      1
# 0       LLNEXT000893    LLNEXT300871    25096550        1.000   0.0000  0.0286  0.0100  0.2353  0.3527  0.0573  0.2184  0.0000  0.1092  3rd     1
# 0       LLNEXT001595    LLNEXT012106    25096550        1.000   0.0000  0.0420  0.0036  0.3956  0.1429  0.2355  0.5098  0.2241  0.4790  FS      1
# 0       LLNEXT004529    LLNEXT006563    25096550        1.000   0.0000  0.0428  0.0039  0.3977  0.1599  0.2328  0.4959  0.1906  0.4385  FS      1
# 0       LLNEXT005203    LLNEXT005302    25096550        1.000   0.0000  0.0429  0.0040  0.3982  0.1652  0.2314  0.4444  0.2072  0.4294  FS      1
# 0       LLNEXT006220    LLNEXT125568    25096550        1.000   0.0000  0.0277  0.0115  0.2276  0.3953  0.0308  0.1524  0.0000  0.0762  4th     0.5
# 0       LLNEXT010135    LLNEXT011783    25096550        1.000   0.0000  0.0282  0.0115  0.2297  0.4031  0.0339  0.1432  0.0000  0.0716  4th     0.5

# Found 6 2nd degree relative pairs (FS)
# Manually selected 6 mothers with the least data layers available: mothers_to_remove_relatives.txt
plink --bfile ../all_chr.mothers --remove mothers_to_remove_relatives.txt --make-bed --out ../all_chr.mothers.no_rel

# PCA
r=0.2
f=all_chr.mothers.no_rel
plink --maf 0.05 --hwe 1e-4 --geno 0.05 --bfile ${f} --out ${f}.flt --make-bed
plink --bfile ${f}.flt --indep-pairwise 1000 50 ${r} --out ${f}.pruned_r2_${r}
plink --bfile $f --extract ${f}.pruned_r2_${r}.prune.in --make-bed --out ${f}.pruned_r2_${r}
plink --pca --bfile ${f}.pruned_r2_${r} --out ${f}.pruned_r2_${r}


#
# Babies
#

# Relatedness IBD analysis

king -b ../all_chr.babies.bed --related --prefix all_chr.babies

grep -v "UN" all_chr.babies.kin
# FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    HetConc HomIBS0 Kinship IBD1Seg IBD2Seg PropIBD InfType Error
# 0       LLNEXT012935    LLNEXT012944    25096550        1.000   0.0000  0.0764  0.0000  0.9963  0.0001  0.4990  0.0067  0.9810  0.9843  Dup/MZ  1
# 0       LLNEXT203282    LLNEXT223335    25096550        1.000   0.0000  0.0750  0.0000  0.9967  0.0001  0.4992  0.0057  0.9830  0.9859  Dup/MZ  1
# 0       LLNEXT204122    LLNEXT205217    25096550        1.000   0.0000  0.0431  0.0034  0.3957  0.1404  0.2382  0.5397  0.1967  0.4666  FS      1
# 0       LLNEXT207937    LLNEXT300081    25096550        1.000   0.0000  0.0462  0.0033  0.4251  0.1352  0.2554  0.5113  0.2226  0.4782  FS      1

# None of these baby ids are in metabolite data, so excluded the first

plink --bfile ../all_chr.babies --remove babies_to_remove_relatives.txt --make-bed --out ../all_chr.babies.no_rel


# PCA
r=0.2
f=all_chr.babies.no_rel
plink --maf 0.05 --hwe 1e-4 --geno 0.05 --bfile ${f} --out ${f}.flt --make-bed 
plink --bfile ${f}.flt --indep-pairwise 1000 50 ${r} --out ${f}.pruned_r2_${r}
plink --bfile $f --extract ${f}.pruned_r2_${r}.prune.in --make-bed --out ${f}.pruned_r2_${r}
plink --pca --bfile ${f}.pruned_r2_${r} --out ${f}.pruned_r2_${r}


# Plot PCs for both mothers and babies
Rscript ${script_dir}/plot_pca.R ${d}


# based on the plots decided to add PC1 and PC2 as covariates for both mothers and babies and to remove 2 PC3 outliers from mothers ()
cat mothers_pca_outliers.txt
# 0	LLNEXT000893
# 0	LLNEXT300871