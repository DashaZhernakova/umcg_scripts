#!/bin/bash
#SBATCH --job-name=SV___BACLIST__
#SBATCH --output=logs/run_GWAS_SV___BACLIST__.out
#SBATCH --error=logs/run_GWAS_SV___BACLIST__.err
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load PLINK

cohort=300OB
type=dSV

d=/data/umcg-tifn/SV/SV_GWAS/
res_dir=${d}/results/${type}/${cohort}/

geno_file=${d}/genotypes/${cohort}/${cohort}_filtered
pheno_file=${d}/data/${cohort}.${type}.filtered.txt
covar_file=${d}/data/${cohort}.covariates.txt

while read line
do
 sv=$line
 echo $sv
 sp=`grep -w "$sv" ${d}/data/sv_name_conversion_table.txt | cut -f4 | uniq`

#
# run real GWAS analysis
#
plink2 \
    --glm sex hide-covar \
    --bfile ${geno_file} \
    --pheno ${pheno_file} \
    --pheno-name "${sv}" \
    --covar ${covar_file} \
    --covar-name "age,read_number,$sp,PC1,PC2" \
    --covar-variance-standardize \
    --out ${res_dir}/${type}.${cohort} \
    --pfilter 0.05

# sort results and convert to EMP
tail -n+2 ${res_dir}/${type}.${cohort}.${sv}.glm.logistic | \
sort -k12,12g | \
python3 ~/scripts/umcg_scripts/SV_GWAS/plink2_assoc_to_EMP.py stdin ${sv} ${cohort} false | tail -n+2 \
> ${res_dir}/${type}.${cohort}.${sv}.eQTLs.txt

# when to compress??

#
# run permuted GWAS
#

done < __BACLIST__

