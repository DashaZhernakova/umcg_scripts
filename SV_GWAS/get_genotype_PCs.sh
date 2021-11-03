module load PLINK

in_f=$1

plink2 \
--bfile $in_f \
--indep-pairwise 1000 50 0.2 \
--out ${in_f}.pruned

plink2 \
--bfile ${in_f} \
--extract ${in_f}.pruned.prune.in \
--make-bed \
--out ${in_f}.pruned

plink2 --bfile ${in_f}.pruned --pca --out ${in_f}.pruned
cut -f1-4  ${in_f}.pruned.eigenvec > ${in_f}.PC1-2.txt
rm ${in_f}.pruned.*