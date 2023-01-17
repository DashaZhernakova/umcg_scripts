timepoints=("M1" "M2" "M3" "M6" "W2")
for t in ${timepoints[@]}
do
 egrep -w "NEXT_ID|${t}" 221125_HMO_data_groups_epitopes_mgml_cleaned_n1272.txt | egrep "first_pregnancy|NEXT_ID" | cut -f2,10- > HMO_${t}.txt
 Rscript ../../../umcg_scripts/misc/transpose.R HMO_${t}.txt HMO_${t}.t.txt
done