#
# Annotate res
#
module load BEDTools

genes=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/resources/ensembl_b37_genes.bed
genes_prot=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS//resources/ensembl_b37_protcoding_genes.bed


timepoints=("M1" "M2" "M3" "M6" "W2")
for t in ${timepoints[@]}
do

  f=/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/plink/HMO_${t}.results.5e-08

  # Add rs ids from dbSNP
  python3 ~/scripts/umcg_scripts/misc/add_rs_by_position.py \
  ${f}.txt /groups/umcg-llnext/tmp01/umcg-dzhernakova/resources/dbsnp_137.b37.vcf.gz 2 \
  > ${f}.rsids.txt


  awk '{FS=OFS="\t"}; {if (NR != 1) print $1, $2-1,$2,$3}' ${f}.rsids.txt > ${f}.tmp.bed

  echo -e "SNPCoord\tOverlapping_genes" > ${f}.tmp.genes
  intersectBed -a ${f}.tmp.bed -b $genes -wa -wb | \
    python2.7 ~/scripts/umcg_scripts/misc/collapseIntersectBedRes.py - \
    >> ${f}.tmp.genes

  python2.7 ~/scripts/umcg_scripts/misc/add_columns_from_file.py  \
    -i ${f}.rsids.txt --header -f ${f}.tmp.genes -i_m 2 -f_m 0 -f_cols 1 \
    > ${f}.tmp.genes1.txt

  echo -e "SNPCoord\tProtcoding_genes_within_250kb" > ${f}.tmp.genes2
  windowBed -a ${f}.tmp.bed -b $genes_prot -w 250000 | \
    python2.7 ~/scripts/umcg_scripts/misc/collapseIntersectBedRes.py - \
    >> ${f}.tmp.genes2

  python2.7 ~/scripts/umcg_scripts/misc/add_columns_from_file.py  \
    -i ${f}.tmp.genes1.txt -f ${f}.tmp.genes2 -i_m 2 -f_m 0 -f_cols 1 --header \
    > ${f%txt}genes.txt


    rm ${f}.tmp*

  python2.7 ~/scripts/umcg_scripts/misc/add_columns_from_file.py \
  -i ${f%txt}genes.txt \
  -f /groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/resources/gwas_catalog_v1.0.2-associations_e105_r2021-12-21.cut.collated.txt \
  -i_m 2 -f_m 0 -f_cols 1 --header \
  > ${f%txt}.genes.gwas_catalog.txt

  rm ${f%txt}genes.txt
  rm *rsids.txt
done
