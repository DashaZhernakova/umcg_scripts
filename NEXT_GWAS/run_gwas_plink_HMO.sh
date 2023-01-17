timepoints=("M1" "M2" "M3" "M6" "W2")
for t in ${timepoints[@]}
do

	plink2 \
		--glm \
		--bfile ../../genotypes/all_chr.mothers.postQC \
		--pheno ../HMO_${t}.INT.txt \
		--covar ../../genotypes/qc/mothers.PC1-2.txt \
		--maf 0.05 --geno 0.05 --hwe 1e-04 \
		--pfilter 0.0001 \
		--out HMO_${t}.results 


	head -1 `ls HMO_${t}.results*linear | head -1` > HMO_${t}.results.5e-08.txt
	
	for f in HMO_${t}.results*linear
	do
		sort -k12,12g $f > ${f}.sorted.txt
		awk '{if ($12 < 5e-08 ) print $0}' ${f}.sorted.txt >> HMO_${t}.results.5e-08.txt
	done

done