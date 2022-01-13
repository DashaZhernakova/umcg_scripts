#!/bin/bash
#SBATCH --job-name=SV___BACLIST__
#SBATCH --output=logs/run_GWAS_SV___BACLIST__.out
#SBATCH --error=logs/run_GWAS_SV___BACLIST__.err
#SBATCH --time=50:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load PLINK
module load Metal

type=dSV
d=/data/umcg-tifn/SV/SV_GWAS/
nperm=10
#sv="B.adolescentis:47"

while read line
do
    sv=$line
    echo $sv
    sp=`grep -w "$sv" ${d}/data/sv_name_conversion_table.txt | cut -f4 | uniq`
    all_nsamples=()

    meta_out_dir=${d}/results/${type}/meta/${sv}/
    meta_out_filebase=${meta_out_dir}/${sv}.meta_res
    mkdir -p ${d}/results/${type}/meta/${sv}/

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

        res_dir=${d}/results/${type}/${cohort}/${sv}/
        mkdir -p ${res_dir}/permutations/

        geno_file=${d}/genotypes/${cohort}/${cohort}_filtered
        perm_fam_dir=${d}/genotypes/${cohort}/permuted/
        pheno_file=${d}/data/${cohort}.${type}.filtered.txt
        covar_file=${d}/data/${cohort}.covariates.txt   
        
        covars="age,read_number,$sp,PC1,PC2"
        if [ $cohort == "DAG3" ]
        then
            echo "DAG3! Use PCs 1-5"
            covars="age,read_number,$sp,PC1,PC2,PC3,PC4,PC5"
        fi

        #
        # run real GWAS analysis
        #
        plink2 \
            --glm sex hide-covar \
            --bfile ${geno_file} \
            --pheno ${pheno_file} \
            --pheno-name "${sv}" \
            --covar ${covar_file} \
            --covar-name "${covars}" \
            --covar-variance-standardize \
            --out ${res_dir}/${type}.${cohort}.${sv} 
        
        echo "$sv real analysis plink return code: $?"
        n=`head -2 ${res_dir}/${type}.${cohort}.${sv}.${sv}.glm.logistic | tail -1 | awk '{print $8}'`
        all_nsamples+=( $n )
        # format the assoc results for METAL: add A1 and A2
        awk '{OFS="\t"}; {if ($6 == $4) {oa=$5}; if ($6 == $5) {oa=$4}; if (NR == 1) {oa="A2"}; {print $1,$2,$3,$6, oa, $7,$8,$9,$10,$11,$12}}' \
        ${res_dir}/${type}.${cohort}.${sv}.${sv}.glm.logistic | gzip -c > ${res_dir}/${type}.${cohort}.${sv}.${sv}.glm.logistic.gz
        reformat_returncode=$?
        
        rm ${res_dir}/${type}.${cohort}.${sv}.${sv}.glm.logistic
        
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
                --covar-name "${covars}" \
                --covar-variance-standardize \
                --out ${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}
                

            echo "$sv permutation $i plink return code: $?"
            # format the assoc results for METAL: add A1 and A2
            awk '{OFS="\t"}; {if ($6 == $4) {oa=$5}; if ($6 == $5) {oa=$4}; if (NR == 1) {oa="A2"}; {print $1,$2,$3,$6, oa, $7,$8,$9,$10,$11,$12}}'  \
            ${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic | gzip -c > ${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic.gz
            reformat_returncode=$?
            
            rm ${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic 
            
            echo -e "PROCESS\t${res_dir}/permutations/${type}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic.gz" >> /data/umcg-tifn/SV/SV_GWAS/scripts/metal_per_sv/${sv}.metal.perm${i}.txt
        done
    done

    #
    # Run meta-analysis
    #
    echo -e "OUTFILE\t${meta_out_filebase} .tbl\nANALYZE\nQUIT" >> $metal_script
    metal $metal_script
    echo "${sv}, real analysis metal return code: $?"

    for i in `seq 1 $nperm`
    do
        echo -e "OUTFILE\t${meta_out_filebase}.perm${i} .tbl\nANALYZE\nQUIT" >>  /data/umcg-tifn/SV/SV_GWAS/scripts/metal_per_sv/${sv}.metal.perm${i}.txt
        metal /data/umcg-tifn/SV/SV_GWAS/scripts/metal_per_sv/${sv}.metal.perm${i}.txt
        echo "${sv}, permutation $i metal return code: $?"
    done
    
    # remove per cohort results
    for cohort in ${cohorts_with_sv[@]}
    do 
        rm ${d}/results/${type}/${cohort}/${sv}/${type}.${cohort}.${sv}.${sv}.glm.logistic.gz
        rm ${d}/results/${type}/${cohort}/${sv}/permutations/${type}.${cohort}.${sv}.perm*.${sv}.glm.logistic.gz
    done

    cohorts_joined=`printf -v var '%s,' "${cohorts_with_sv[@]}"; echo "${var%,}"`
    samplesize_joined=`printf -v var '%s,' "${all_nsamples[@]}"; echo "${var%,}"`


    # sort results and convert to EMP per sv!
    tail -n+2 ${meta_out_filebase}1.tbl | \
    sort -k6g | \
    python3 ~/scripts/umcg_scripts/SV_GWAS/metal_to_EMP.py stdin ${sv} $cohorts_joined $samplesize_joined 0.05 | tail -n+2 \
    > ${meta_out_filebase}.eQTLs.txt
    
    for p in `seq 1 $nperm`
    do
        tail -n+2 ${meta_out_filebase}.perm${p}1.tbl | \
        sort -k6g | \
        python3 ~/scripts/umcg_scripts/SV_GWAS/metal_to_EMP.py stdin ${sv} $cohorts_joined $samplesize_joined 0.05 | tail -n+2 \
        > ${meta_out_filebase}.eQTLs.perm${p}.txt
        rm ${meta_out_filebase}.perm${p}1.tbl
    done
    rm ${meta_out_filebase}1.tbl

done < __BACLIST__