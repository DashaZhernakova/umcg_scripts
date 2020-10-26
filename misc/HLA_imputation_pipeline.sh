#
### Alignment of phased chr6 genotypes to the HLA reference
#

## Uses Genotype Harmonizer: https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer-Download

## variables to edit
# Path to chr 6 genotypes phased by SHAPEIT2:
input=5_phaseing/chr_6
# Path to the output aligned chr 6 genotypes:
output=X5_alignPhasingToHLA/chr_6
# HLA reference panel in PLINK format:
HLAReferencePlink=/groups/umcg-wijmenga/tmp04/umcg-mjbonder/HLAImputeReference/B37/HLA_T1DGC/T1DGC_REF_B37

java -Xms50g -Xmx50g -jar GenotypeHarmonizer-1.4.9/GenotypeHarmonizer.jar \
    -i ${input} \
    -I SHAPEIT2 \
    -o ${output} \
    -O SHAPEIT2 \
    -r ${HLAReferencePlink} \
    -R PLINK_BED \
    --update-id \
    --variants 500 \
    --mafAlign 0.01 \
    --min-ld 0.4 \
    --forceChr 6 \
    --min-variants 3


#
### Imputation
#

## variables to edit

# aligned phased haplotypes from the previous step:
known_haps_g=${output}.haps

#Path to the resulting output file
output_imputed="./imputed/chr6.HLA_imputed"

# genetic map downloaded from https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html :
m="/gcc/resources/b37/imputation/1000G/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr6_combined_b37.txt"

# HLA reference panel in impute2 format:
h="/groups/umcg-wijmenga/tmp04/umcg-mjbonder/HLAImputeReference/B37/HLA_T1DGC/Imput_T1DGC_REF_B37.hap.gz"
l="/groups/umcg-wijmenga/tmp04/umcg-mjbonder/HLAImputeReference/B37/HLA_T1DGC/Imput_T1DGC_REF_b37.legend.gz"

# Additional impute2 parameters
additonalImpute2Param="-Ne 20000 -k_hap 12000"

impute2 \
        -known_haps_g $known_haps_g \
        -m $m \
        -h $h \
        -l $l \
        -o $output_imputed \
        -use_prephased_g \
        $additonalImpute2Param


