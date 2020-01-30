#!/usr/local/bin/python

'''
Adds rs ids and reformats the tables for MR

Usage example:
#python /groups/umcg-lld/tmp03/umcg-dzhernakova/MR/get_rs_by_positions_dbSNP_adriaan_MR.py \
#/groups/umcg-wijmenga/tmp03/users/umcg-avdgraaf/ll_phewas/2019_09_23_analysis/beta_se_files/ \
#/groups/umcg-wijmenga/tmp03/users/umcg-avdgraaf/ll_phewas/2019_09_23_analysis/filtered/ \
#__BACT__ \
#/groups/umcg-lld/tmp03/umcg-dzhernakova/All_20180423.vcf.gz \
0  \
> /groups/umcg-lld/tmp03/umcg-dzhernakova/MR/data/mibiogenSep2019/__BACT__.filtered.rsids.txt
'''




import sys
import os
import vcf 
import gzip
from collections import defaultdict

def getRs(chr_pos, vcf_reader):
    chr = chr_pos.split(":")[0]
    pos = int(chr_pos.split(":")[1])
    for record in vcf_reader.fetch(chr, pos - 1, pos):
        #print "##", record.ID, record
        if record.is_snp:
            return record.ID
    return chr_pos

beta_file_dir = sys.argv[1]
orig_file_dir = sys.argv[2]
pheno = sys.argv[3]
vcf_reader = vcf.Reader(open(sys.argv[4]))
col_num = int(sys.argv[5])

snp2rs = defaultdict()

print "SNP\tbeta\tse\tsamplesize\teffect_allele\tother_allele\teaf\tPhenotype\tpval"
in_fname = beta_file_dir + "/" + pheno + ".filtered.txt"
orig_f = gzip.open(orig_file_dir + "/" + pheno + ".filtered.txt.gz")
with open(in_fname) as f:
    f.readline()
    orig_f.readline()
    for l in f:
        spl = l.strip().split("\t")
	spl_orig = orig_f.readline().strip().split("\t")
	if not spl_orig[1] == spl[col_num]:
		print "SNPs differ:", spl_orig[1], spl[0]
		sys.exit(1)
       	snp = spl[col_num]
        rs = snp2rs.get(snp, None)
        if not rs:
     		rs = getRs(snp, vcf_reader)
            	snp2rs[snp] = rs
        spl[col_num] = rs
	pval = spl_orig[0]
	if spl[6] == "nan":
		spl[6] = "0.5"
        print "\t".join(spl) + "\t" + pheno + "\t" + pval
