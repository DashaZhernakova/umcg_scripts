#!/bin/bash
#SBATCH --job-name=fdr
#SBATCH --output=fdr.out
#SBATCH --error=fdr.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=85gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

# merge sorted files
meta_comb_dir=/data/umcg-tifn/SV/SV_GWAS/results/dSV/meta_combined/
nperm=10

mkdir ${meta_comb_dir}
echo -e "PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR" \
> ${TMPDIR}/eQTLs.txt

#
sort -m -k1,1g -T $TMPDIR /data/umcg-tifn/SV/SV_GWAS/results/dSV/meta/*/*.meta_res.eQTLs.txt \
>> ${TMPDIR}/eQTLs.txt

head -1000001 ${TMPDIR}/eQTLs.txt | gzip -c > ${meta_comb_dir}/eQTLs.txt.gz


#meta_comb_dir=/data/umcg-tifn/SV/SV_GWAS/results/dSV/meta_combined/

for p in `seq 1 $nperm`
do
    echo -e "PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR" \
    > ${TMPDIR}/PermutedEQTLsPermutationRound${p}.txt

    sort -m -k1,1g -T $TMPDIR /data/umcg-tifn/SV/SV_GWAS/results/dSV/meta/*/*.meta_res.eQTLs.perm${p}.txt \
    >> ${TMPDIR}/PermutedEQTLsPermutationRound${p}.txt

    head -1000001 ${TMPDIR}/PermutedEQTLsPermutationRound${p}.txt | gzip -c > ${meta_comb_dir}/PermutedEQTLsPermutationRound${p}.txt.gz
done


java -Xmx80g -Xms80g -jar /data/umcg-tifn/SV/SV_GWAS/eqtl-mapping-pipeline-1.3.9-SNAPSHOT/eqtl-mapping-pipeline.jar \
--mode util \
--fdr \
--in  ${meta_comb_dir} \
--threshold 0.05 \
--perm 10 --nreqtls 1000000