module load plink

ts=("MTL_CD45+_CD20-" "MTL_CD45-" "MTL_gran" "MTL_CD57+" "MTL_CD20+" "MTL_lymp")
cutoff=0.05

for t in ${ts[@]}
do
  echo $t
  zcat ../eQTLs.txt.gz | head -1 > ${t}.QTLs.p${cutoff}.txt
  zgrep $t ../eQTLs.txt.gz | awk -v p="$cutoff" 'BEGIN {FS=OFS="\t"}; {if ($1 < p) print $0}' >> ${t}.QTLs.p${cutoff}.txt


  plink \
    --bfile /groups/umcg-lld/scr01/dasha/genotypes/combined_plink \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump ${t}.QTLs.p${cutoff}.txt \
    --clump-snp-field SNPName \
    --clump-field PValue \
    --out ${t}.QTLs.p${cutoff}

  awk 'NR!=1{print $3}' ${t}.QTLs.p${cutoff}.clumped >  ${t}.QTLs.p${cutoff}.clumped.snp

  #awk '{print $1,$8}' ${t}.QTLs.p${cutoff}.txt > ${t}.QTLs.p${cutoff}.SNP.pvalue

  plink \
    --bfile /groups/umcg-lld/scr01/dasha/genotypes/combined_plink \
    --score ${t}.QTLs.p${cutoff}.txt 2 10 18 header \
    --extract ${t}.QTLs.p${cutoff}.clumped.snp \
    --out ${t}.QTLs.p${cutoff}

done
rm *nosex
rm *clumped*
rm *nopred
gzip *.QTLs.p${cutoff}.txt
