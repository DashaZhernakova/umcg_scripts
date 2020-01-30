#!/usr/bin/python
import sys
import gzip
import re

maf_dict = {}
snp_f = gzip.open(sys.argv[2])
snp_f.readline()
snp_f.readline()
for l in snp_f:
	if not l.startswith("SNP alleles"):
		spl = l.strip().split("\t")
		#print spl
		maf_dict[spl[0]] = spl[7]
snp_f.close()

sys.stderr.write("Read MAFs: " + str(len(maf_dict)) + "\n")

qtls_type = "NA"
#if len(sys.argv) > 3:
#	qtls_type = sys.argv[3]
if sys.argv[1].endswith(".gz"):
	f = gzip.open(sys.argv[1])
else:
	f = open(sys.argv[1])
f.readline()
print "Phenotype\tSNP\tbeta\tse\tpval\teffect_allele\tother_allele\ttype\teaf\tchr:pos\tsamplesize"
for l in f:
	spl = l.strip().split("\t")
	#print spl
	beta = spl[19]
	if ";" in beta:
		if beta.startswith("-;"):
			m = re.match("-;(.*)\s+\((.*)\)", beta)
		else:
			m = re.match("(.*)\s+\((.*)\);", beta)
	else:
		m = re.match("(.*)\s+\((.*)\)", beta)	
	
	a2 = set(spl[8].split("/")) - set(spl[9])
	chrpos = spl[2] + ":" + spl[3]
	maf = maf_dict[spl[1]]
	qtls_type = spl[7]
	#print spl[1], spl[16], spl[19]
	#print m
	#print m.group(0)
	#print m.group(1)
	#print m.group(2)
	print spl[16].split("_")[1].replace(".", "-") + "\t" + spl[1] + "\t" + m.group(1) + "\t" + m.group(2) + "\t" + spl[0] + "\t" + spl[9] + "\t" + a2.pop() + "\t" + qtls_type + "\t" + maf + "\t" + chrpos + "\t1178" 
