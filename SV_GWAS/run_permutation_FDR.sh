

sv=Blautia__wexlerae__DSM__19850:236
mkdir EMP
python metal_processed_to_EMP.py ../meta/dSV_meta_all_1e-06.sorted.rsids.txt | gzip -c > EMP/eQTLs.txt.gz

for p in `seq 1 100`
do
   python metal_to_EMP.py meta/meta_perm_${p}_1.tbl ${sv} LLD,300OB \
   | (read h; echo "$h"; sort -k1,1g) \
   | gzip -c > EMP/PermutedEQTLsPermutationRound${p}.txt.gz
done

java -Xmx45g -Xms45g -jar /groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/tools/eqtl-mapping-pipeline-1.3.9-SNAPSHOT/eqtl-mapping-pipeline.jar \
--mode util --fdr \
--in EMP/ \
--threshold 0.05 \
--perm 100 \
--nreqtls 100000


cutoff=1e-04
res_f="${base_dir}/results/meta/dSV_meta_all_${cutoff}.txt"
echo -e "SV\tMarkerName\tAllele1\tAllele2\tWeight\tZscore\tP-value\tDirection\tCohort" > $res_f
while read line
do
 sv=$line
 echo $sv
 cohort="LLD,300OB"
 zcat ../results/meta/dSV_meta_${sv}_1.tbl.gz | awk -v bac=${sv} -v d=${cohort} -v p=${cutoff} '{OFS="\t"}; {if ($6 < p) {print bac, $0, d}}'  >> $res_f
done < all_bacs_LLD+300OB.txt