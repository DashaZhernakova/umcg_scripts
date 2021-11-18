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
module load Metal
set -e

type=dSV
d=/data/umcg-tifn/SV/SV_GWAS/
cohorts_nodag3=("300OB" "500FG" "LLD")
nperm=5

while read line
do
    sv=$line
    echo $sv
    sp=`grep -w "$sv" ${d}/data/sv_name_conversion_table.txt | cut -f4 | uniq`
    all_nsamples=()

    meta_out_filebase=${d}/results/${type}/meta/${sv}.meta_res
    
    metal_script=/data/umcg-tifn/SV/SV_GWAS/scripts/metal_per_sv/${sv}.metal.txt
    cat /data/umcg-tifn/SV/SV_GWAS/scripts/metal_header.txt > $metal_script
    for p in `seq 1 $nperm`
    do
       cat /data/umcg-tifn/SV/SV_GWAS/scripts/metal_header.txt > /data/umcg-tifn/SV/SV_GWAS/scripts/metal_per_sv/${sv}.metal.perm${p}.txt
    done

    # check in which cohorts the SV is present
    IFS=',' read -ra cohorts_with_sv <<< `grep -w $sv ${d}/data/dSV_per_cohort.txt | cut -f6`

    for cohort in ${cohorts_with_sv[@]}
    do
        echo -e "\n\nRunning the analysis for ${cohort}\n\n"

        res_dir=${d}/results/${type}/${cohort}/
        mkdir -p ${res_dir}/permutations/

        geno_file=${d}/genotypes/${cohort}/${cohort}_filtered
        perm_fam_dir=${d}/genotypes/${cohort}/permuted/
        pheno_file=${d}/data/${cohort}.${type}.filtered.txt
        covar_file=${d}/data/${cohort}.covariates.txt   
        
        
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
            --out ${res_dir}/${type}.${cohort}.${sv}
        
        plink_returncode=$?
        n=`head -2 ${res_dir}/${type}.${cohort}.${sv}.${sv}.glm.logistic | tail -1 | awk '{print $7}'`
        all_nsamples+=( $n )
        # format the assoc results for METAL: add A1 and A2
        awk '{OFS="\t"}; {if ($6 == $4) {oa=$5}; if ($6 == $5) {oa=$4}; if (NR == 1) {oa="A2"}; {print $1,$2,$3,$6, oa, $7,$8,$9,$10,$11,$12}}' \
        ${res_dir}/${type}.${cohort}.${sv}.${sv}.glm.logistic | gzip -c > ${res_dir}/${type}.${cohort}.${sv}.${sv}.glm.logistic.gz
        reformat_returncode=$?
        
        if [ $plink_returncode -eq 0 ] && [ $reformat_returncode -eq 0 ]
        then 
            rm ${res_dir}/${type}.${cohort}.${sv}.${sv}.glm.logistic
        else
            echo "ERROR in GWAS or awk reformatting. plink_returncode=$plink_returncode ; reformat_returncode=$reformat_returncode . Output file: ${res_dir}/${type}.${cohort}.${sv}.${sv}.glm.logistic "
        fi
        echo -e "PROCESS\t${res_dir}/${type}.${cohort}.${sv}.${sv}.glm.logistic.gz\n" >> $metal_script
        
        
        
        #
        # run permuted GWAS
        #
        for i in `seq 1 $nperm`
        do
            plink2 \
                --glm sex hide-covar \
                --bed ${geno_file}.bed \
                --bim ${geno_file}.bim \
                --fam ${perm_fam_dir}/perm${i}.fam \
                --pheno ${pheno_file} \
                --pheno-name "${sv}" \
                --covar ${covar_file} \
                --covar-name "age,read_number,$sp,PC1,PC2" \
                --covar-variance-standardize \
                --out ${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}
            plink_returncode=$?

            # format the assoc results for METAL: add A1 and A2
            awk '{OFS="\t"}; {if ($6 == $4) {oa=$5}; if ($6 == $5) {oa=$4}; if (NR == 1) {oa="A2"}; {print $1,$2,$3,$6, oa, $7,$8,$9,$10,$11,$12}}'  \
            ${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic | gzip -c > ${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic.gz
            reformat_returncode=$?
            
            if [ $plink_returncode -eq 0 ] && [ $reformat_returncode -eq 0 ]
            then 
                rm ${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic 
            else
                echo "ERROR in GWAS or awk reformatting. plink_returncode=$plink_returncode ; reformat_returncode=$reformat_returncode . Output file: ${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic  "
            fi  
            
            echo -e "PROCESS\t${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic.gz" >> /data/umcg-tifn/SV/SV_GWAS/scripts/metal_per_sv/${sv}.metal.perm${i}.txt
        done
    done

    #
    # Run meta-analysis
    #
    echo -e "OUTFILE\t${meta_out_filebase} .tbl\nANALYZE\nQUIT" >> $metal_script
    metal $metal_script

    for i in `seq 1 $nperm`
    do
        #echo -e "OUTFILE\t${meta_out_filebase}.perm${i} .tbl\nANALYZE\nQUIT" >>  /data/umcg-tifn/SV/SV_GWAS/scripts/metal_per_sv/${sv}.metal.perm${i}.txt
        metal /data/umcg-tifn/SV/SV_GWAS/scripts/metal_per_sv/${sv}.metal.perm${i}.txt
    done
    
done < __BACLIST__

B.longum:17
E.rectale:66
B.wexlerae:70
D.longicatena:150
D.formicigenerans:364

cohorts_joined=`printf -v var '%s,' "${cohorts_with_sv[@]}"; echo "${var%,}"`
samplesize_joined=`printf -v var '%s,' "${all_nsamples[@]}"; echo "${var%,}"`


#PermutedEQTLsPermutationRound10.txt.gz
# sort results and convert to EMP
tail -n+2 ${meta_out_filebase}1.tbl | \
sort -k6g | \
python3 ~/scripts/umcg_scripts/SV_GWAS/metal_to_EMP.py stdin ${sv} $cohorts_joined $samplesize_joined | tail -n+2 \
> ${d}/results/${type}/meta/${sv}.meta_res.eQTLs.txt
#rm ${res_dir}/${type}.${cohort}.${sv}.glm.logistic