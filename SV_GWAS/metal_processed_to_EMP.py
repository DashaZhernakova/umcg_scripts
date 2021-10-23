import sys
import gzip

empty_out_spl=["1","x","1","1","x","1","1","trans","A/T","A","0","x","0","0","0","0","x","0","0","0","0","0"]

fname = sys.argv[1]
max_lines = 100000

if fname.endswith(".gz"):
    f = gzip.open(fname)
else:
    f = open(fname)

f.readline()
print("PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR")

line_num = 1

for l in f:
    spl = l.strip().split()
    snp = spl[2]
    p = float(spl[7])
    #if p > 0.05:
    #    continue
    if line_num > max_lines:
        break
    zscore = spl[6]
    empty_spl_cp = empty_out_spl[:]
    empty_spl_cp[0] = str(p)
    empty_spl_cp[1] = spl[1]
    empty_spl_cp[2] = snp.split(":")[0]
    empty_spl_cp[3] = snp.split(":")[1]
    empty_spl_cp[4] = spl[0]
    empty_spl_cp[8] = spl[3].upper() + "/" + spl[4].upper()
    empty_spl_cp[9] = spl[3].upper()
    empty_spl_cp[16] = spl[0] 
    empty_spl_cp[10] = zscore
    empty_spl_cp[11] = spl[9]
    empty_spl_cp[12] = zscore
    empty_spl_cp[13] = spl[5]

    print ("\t".join(empty_spl_cp))
    line_num += 1





