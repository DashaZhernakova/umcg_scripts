ml PLINK/1.9-beta6-20190617

for d in /groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/EMP_HMO_*/
do
    cd $d
    echo "Starting with $d"
    hmos=( $( cut -f5 eQTLs.1e-5.txt | tail -n+2 | sort | uniq ) )
    mkdir clumped
    for hmo in ${hmos[@]}
    do
     awk -v h=$hmo '{OFS=FS="\t"}; {if ($5 == h || NR == 1) print}' eQTLs.1e-5.txt > clumped/tmp.${hmo}.eQTLs.1e-5.txt
     plink --bfile /groups/umcg-llnext/tmp01/umcg-dzhernakova/genotypes/1000G_snps_nodup \
     --clump clumped/tmp.${hmo}.eQTLs.1e-5.txt \
     --clump-snp-field SNPName \
     --clump-field PValue \
     --out clumped/${hmo}.1e-5 \
     --clump-r2 0.4
     
     awk -v h="${hmo}" '{if ($3 != "") print $3 "\t" h}' clumped/${hmo}.1e-5.clumped | tail -n+2 > clumped/${hmo}.1e-5.clumped.SNP-HMO.txt
    done
    cd /groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/
done