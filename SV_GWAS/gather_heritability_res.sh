while read line
do
	sv=$line
	h=`grep "V(G)/Vp" ../results/vSV/heritability/${sv}/${sv}_norel.hsq | awk '{print $2,$3}'`
	p=`grep "Pval" ../results/vSV/heritability/${sv}/${sv}_norel.hsq | awk '{print $2}'`
	echo "$sv $h $p"
done < dag3_vSVs.txt



while read line
do
	sv=$line
	h=`grep "Sum of V(G)/Vp" ../results/vSV/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $4,$5}'`
	p=`grep "Pval" ../results/vSV/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	echo "$sv $h $p"
done < dag3_vSVs.txt