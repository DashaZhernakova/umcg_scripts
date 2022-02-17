cd /data/umcg-tifn/SV/SV_GWAS/genotypes
snp="9:136149830"
echo $snp > snps_for_extraction/${snp}.txt

cohorts=("LLD" "500FG" "DAG3")
for c in ${cohorts[@]}
do
    mkdir ${c}/text_genotypes
    java -Xmx10g -jar ~/tools/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar \
    -i ${c}/${c}_filtered -I PLINK_BED \
    -O table -o  ${c}/text_genotypes/${c}.${snp} \
    -vf snps_for_extraction/${snp}.txt
done