ml PLINK
timepoints=("Mother_B" "Mother_P12" "Mother_P28")
for t in ${timepoints[@]}
do

	plink2 \
		--glm \
		--bfile ../../genotypes/all_chr.mothers.postQC \
		--pheno metabo_${t}.filtered.txt \
		--covar ../../genotypes/qc/mothers.PC1-2.txt \
		--maf 0.05 --geno 0.05 --hwe 1e-04 \
		--pfilter 0.0001 \
		--out metabo_${t}.plink 


	head -1 `ls metabo_${t}.plink*linear | head -1` > metabo_${t}.plink.5e-08.txt
	
	for f in metabo_${t}.plink*linear
	do
		sort -k12,12g $f > ${f}.sorted.txt
		awk '{if ($12 < 5e-08 ) print $0}' ${f}.sorted.txt >> metabo_${t}.plink.5e-08.txt
	done

done


# Babies

t="Baby_B"

plink2 \
	--glm sex \
	--bfile ../../genotypes/all_chr.babies.postQC \
	--pheno metabo_Baby_B.filtered.txt \
	--covar ../../genotypes/qc/babies.PC1-2.txt \
	--maf 0.05 --geno 0.05 --hwe 1e-04 \
	--pfilter 0.0001 \
	--out metabo_${t}.plink 


head -1 `ls metabo_${t}.plink*linear | head -1` > metabo_${t}.plink.5e-08.txt

for f in metabo_${t}.plink*linear
do
	sort -k12,12g $f > ${f}.sorted.txt
	awk '{if ($12 < 5e-08 ) print $0}' ${f}.sorted.txt >> metabo_${t}.plink.5e-08.txt
done