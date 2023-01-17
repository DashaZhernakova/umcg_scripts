#
# Annotate res
#
module load BEDTools
f=EMP_HMO_M3.eQTLs.5e-08

genes=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/resources/ensembl_b37_genes.bed
genes_prot=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS//resources/ensembl_b37_protcoding_genes.bed


# Extract all associations with P < 5e-08
zcat eQTLs.txt.gz | awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1 || $1 < 5e-8) print }' | cut -f1-5,9-12,14 > ${f}.txt

# Add rs ids from dbSNP
python3 ~/scripts/umcg_scripts/SV_GWAS/clean/gwas_scripts_misc/add_rs_by_position.py \
${f}.txt /groups/umcg-llnext/tmp01/umcg-dzhernakova/resources/dbsnp_137.b37.vcf.gz 1 \
> ${f}.rsids.txt


awk '{FS=OFS="\t"}; {if (NR != 1) print $4, $5-1,$5,$2}' ${f}.rsids.txt > ${f}.tmp.bed

echo -e "SNPCoord\tOverlapping_genes" > ${f}.tmp.genes
intersectBed -a ${f}.tmp.bed -b $genes -wa -wb | \
  python2.7 ~/scripts/umcg_scripts/misc/collapseIntersectBedRes.py - \
  >> ${f}.tmp.genes

 python2.7 ~/scripts/umcg_scripts/misc/add_columns_from_file.py  \
  -i ${f}.rsids.txt --header -f ${f}.tmp.genes -i_m 1 -f_m 0 -f_cols 1 \
  > ${f}.tmp.genes1.txt

echo -e "SNPCoord\tProtcoding_genes_within_250kb" > ${f}.tmp.genes2
windowBed -a ${f}.tmp.bed -b $genes_prot -w 250000 | \
  python2.7 ~/scripts/umcg_scripts/misc/collapseIntersectBedRes.py - \
  >> ${f}.tmp.genes2

 python2.7 ~/scripts/umcg_scripts/misc/add_columns_from_file.py  \
  -i ${f}.tmp.genes1.txt -f ${f}.tmp.genes2 -i_m 1 -f_m 0 -f_cols 1 --header \
  > ${f%txt}genes.txt


  rm ${f}.tmp*

python2.7 ~/scripts/umcg_scripts/misc/add_columns_from_file.py \
-i ${f%txt}genes.txt \
-f /groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/resources/gwas_catalog_v1.0.2-associations_e105_r2021-12-21.cut.collated.txt \
-i_m 1 -f_m 0 -f_cols 1 --header \
> ${f%txt}genes.gwas_catalog.txt


