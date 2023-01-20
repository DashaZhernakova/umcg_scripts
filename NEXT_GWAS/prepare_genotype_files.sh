module load BCFtools
module load PLINK

script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
d=/groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/
# Remove SNPs with imputation quality < 0.5
for chr in `seq 1 22`
do
 plink2 \
 --vcf /groups/umcg-llnext/tmp01/data/imputed/rawdata/LLnextv2_2920merged23122022.vcfs/${chr}.pbwt_reference_impute.vcf.gz \
 --extract-if-info "INFO > 0.5" \
 --make-bed --out ${d}/${chr}.info_filtered
 
 awk 'BEGIN {FS=OFS="\t"}; {if ($2 == ".") $2 = $1 ":" $4; print $0}' ${d}/${chr}.info_filtered.bim > ${d}/${chr}.info_filtered.bim.tmp
 mv ${d}/${chr}.info_filtered.bim.tmp ${d}/${chr}.info_filtered.bim
 
 awk 'BEGIN {FS=OFS="\t"}; {$2 = $1 ":" $4; print $0}' ${d}/${chr}.info_filtered.bim > ${d}/${chr}.info_filtered.bim.tmp
 mv ${d}/${chr}.info_filtered.bim.tmp ${d}/${chr}.info_filtered.bim
 
done


# Merge all chromosomes
rm ${d}/merge_list.txt
for chr in `seq 1 22`
do
 echo "${d}/${chr}.info_filtered" >> ${d}/merge_list.txt
done


module load PLINK/1.9-beta6-20190617
plink --merge-list ${d}/merge_list.txt --make-bed --out ${d}/all_chr.info_filtered

cat ${d}/all_chr.filtered-merge.missnp > ${d}/multiallelic_SNPs.txt

for chr in `seq 1 22`
do
 plink --bfile ${d}/${chr}.info_filtered \
 --exclude ${d}/multiallelic_SNPs.txt \
 --make-bed --out ${d}/${chr}.info_filtered.nodup
done

rm ${d}/merge_list.txt
for chr in `seq 1 22`
do
 echo "${d}/${chr}.info_filtered.nodup" >> ${d}/merge_list.txt
done

plink --merge-list ${d}/merge_list.txt --make-bed --out ${d}/intermediate_files/all_chr.info_filtered.all_samples

# Keep only samples with LLNEXT id
plink --bfile ${d}/intermediate_files/all_chr.info_filtered.all_samples --keep ${d}/llnext_samples.txt --make-bed --out ${d}/all_chr.info_filtered

# Separate mothers and babies
plink --bfile ${d}/all_chr.info_filtered.all_samples --keep ${d}/qc/mothers.txt --make-bed --out ${d}/all_chr.mothers
plink --bfile ${d}/all_chr.info_filtered.all_samples --keep ${d}/qc/babies.txt --make-bed --out ${d}/all_chr.babies

# Run genotype QC
bash ${script_dir}/NEXT_genotype_QC.sh
# removed relatives and for mothers PC3 outliers
# filtered genotypes written to *postQC files

# add sex for babies
plink2 --bfile ${d}/all_chr.babies.postQC --update-sex ${d}/baby_sex.txt --make-bed --out ${d}/all_chr.babies.postQC.2
mv ${d}/all_chr.babies.postQC.2.fam ${d}/all_chr.babies.postQC.fam
mv ${d}/all_chr.babies.postQC.2.bim ${d}/all_chr.babies.postQC.bim
mv ${d}/all_chr.babies.postQC.2.bed ${d}/all_chr.babies.postQC.bed


# Convert to Trityper
#java  -Xms50g -Xmx50g -jar ~/tools/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar -i ${d}/all_chr.mothers.postQC -I PLINK_BED -o ${d}/NEXT_mothers -O TRITYPER
#java  -Xms50g -Xmx50g -jar ~/tools/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar -i ${d}/all_chr.babies.postQC -I PLINK_BED -o ${d}/NEXT_babies -O TRITYPER