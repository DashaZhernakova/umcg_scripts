import sys
import scipy.stats as st
import math
import gzip

empty_out_spl=["1","x","1","1","x","1","1","trans","A/T","A","0","x","0","0","0","0","x","0","0","0","0","0"]

fname = sys.argv[1]
pheno = sys.argv[2]
dataset = sys.argv[3]
if len(sys.argv > 4):
    read_header = sys.argv[4]
else:
    read_header = True
#max_lines = 1000000

if fname.endswith(".gz"):
    f = gzip.open(fname)
elif (fname == "stdin"):
    f = sys.stdin
else:
    f = open(fname)
if read_header.lower() == "true":
    f.readline()

print("PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR")

line_num = 1

for l in f:
    spl = l.strip().split()
    p = float(spl[10])
    #if p > 0.05:
    #    continue
    #if line_num > max_lines:
    #    break
    zscore = spl[9]
    empty_spl_cp = empty_out_spl[:]
    empty_spl_cp[0] = str(p)
    empty_spl_cp[1] = spl[2]
    empty_spl_cp[2] = spl[0]
    empty_spl_cp[3] = spl[1]
    empty_spl_cp[4] = pheno
    empty_spl_cp[8] = spl[3] + "/" + spl[4]
    empty_spl_cp[9] = spl[3]
    empty_spl_cp[16] = pheno 
    empty_spl_cp[10] = zscore
    empty_spl_cp[11] = dataset
    empty_spl_cp[12] = zscore
    empty_spl_cp[13] = spl[6]

    print ("\t".join(empty_spl_cp))
    line_num += 1





