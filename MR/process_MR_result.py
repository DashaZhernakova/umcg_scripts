#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import math
from collections import defaultdict
from functools import cmp_to_key
import itertools
import statsmodels.stats.multitest as multi

def isNaN(string):
    return string != string

taxa_to_remove = ["family.Bifidobacteriaceae.id.433", "family.unknownfamily.id.1000001214", "genus.unknowngenus.id.1000001215"]
fname = r'/Users/Dasha/work/UMCG/data/MR/results2/mibiogen/mibiogenOct2020/MR_mibiogen_oct2020.xlsx'
sh1 = 'gwas-mb'
sh2 = 'mb-gwas'
#sh1 = "UKB-mb"
#sh2 = "mb-UKB"
df1 = pd.read_excel (fname, sheet_name=sh1)
df2 = pd.read_excel (fname, sheet_name=sh2)
df1['exposure'] = [x for (x,y) in df1['exposure'].str.split(" \|\| ").to_list()]
df2['outcome'] = [x for (x,y) in df2['outcome'].str.split(" \|\| ").to_list()]

df1['tmp'] = df1['id.exposure'] + ":" + df1['outcome']
df2['tmp'] = df2['id.outcome'] + ":" + df2['exposure']

df1new = pd.merge(df1, df2[['tmp', 'pval']], on ='tmp', how ='left')
df1new.rename(columns = {"pval_x" : "pval", "pval_y" : "reverse_MR_pval"}, inplace=True)

df2new = pd.merge(df2, df1[['tmp', 'pval']], on ='tmp', how ='left')
df2new.rename(columns = {"pval_x" : "pval", "pval_y" : "reverse_MR_pval"}, inplace=True)

df1new.drop('tmp', inplace=True, axis=1)
df2new.drop('tmp', inplace=True, axis=1)

# Write all filters that failed into one field
def get_failed_flts(row):
    flt_line = "-"
    if int(row['nsnp']) < 3:
        flt_line += ';Number of SNPs'
    if not isNaN(row['egger_intercept_pval']):
        if float(row['egger_intercept_pval']) < 0.05:
            flt_line += ';Egger intercept p-value'
    if not isNaN(row['weighted_median_pval']):        
        if float(row['weighted_median_pval']) > 0.05:
            flt_line += ';Weighted median p-value'
    if not isNaN(row['mr_presso_global']):        
        if not isNaN(row['mr_presso_outlier_cor_pval']):
            if row['mr_presso_global'] == '<1e-04':
                if float(row['mr_presso_outlier_cor_pval']) > 0.05:
                    flt_line += ';MR PRESSO outlier test'
            elif float(row['mr_presso_global']) < 0.05 and float(row['mr_presso_outlier_cor_pval']) > 0.05:
                flt_line += ';MR PRESSO outlier test'
    if not isNaN(row['mr_presso_pval']):
        if float(row['mr_presso_pval']) > 0.05:
            flt_line += ';MR PRESSO p-value'
    # Leave-one-out updated to filter out only cases with exactly 1 p-value > 0.05
    loo_pvals = row['leave_one_out_pval']
    pval_cnt = 0
    if not isNaN(loo_pvals):
        for pval in map(float, loo_pvals.split(",")):
            if pval > 0.05:
                pval_cnt += 1
        if pval_cnt == 1:
            flt_line += ';Leave-one-out analysis'
        if not isNaN(row['reverse_MR_pval']):
            if float(row['reverse_MR_pval']) < 0.05:
                flt_line += ';Reverse MR p-value'
    return(flt_line.replace("-;", "", 1))


df1new['failed_filters'] = df1new.apply (lambda row: get_failed_flts(row), axis=1)
df2new['failed_filters'] = df2new.apply (lambda row: get_failed_flts(row), axis=1)


## Filter by samplesize: take the largest of the duplicated traits
samplesize_dict = {}
with open("C:/Users/Dasha/work/UMCG/data/MR/data/MRBase_all_outcomes.271020.txt") as f:
    f.readline()
    for l in f:
        spl = l.split("\t")
        samplesize_dict[spl[1]] = spl[14]

df1new['samplesize'] = df1new['id.exposure'].replace(samplesize_dict, inplace = False)
df2new['samplesize'] = df2new['id.outcome'].replace(samplesize_dict, inplace = False)


df1new = add_samplesize_filter(df1new)
df2new = add_samplesize_filter(df2new)

# cmp 2 MR results first by sample size, then by number of shared SNPs
def comparator(x, y):
    #ss_col = 3
    #nsnp_col = 5
    ss_col = 'samplesize'
    nsnp_col = 'nsnp'
    if int(x[ss_col]) > int(y[ss_col]):
        return -1
    if int(x[ss_col]) < int(y[ss_col]):
        return 1
    if int(x[ss_col]) == int(y[ss_col]):
        if int(x[nsnp_col]) > int(y[nsnp_col]):
            return -1
        if int(x[nsnp_col]) < int(y[nsnp_col]):
            return 1
        if int(x[nsnp_col]) == int(y[nsnp_col]):
            return 0


def add_samplesize_filter(df):
    df['filter_samplesize'] = "FALSE"
    pair2line = defaultdict(list)
    for index, row in df.iterrows():
        pair2line[row['exposure'] + ":" + row['outcome']].append(row)
    for p, row_lst in pair2line.items():
        row_lst.sort(key=cmp_to_key(comparator))
        row_lst[0]['filter_samplesize'] = "TRUE"
    return (pd.DataFrame(list(itertools.chain(*pair2line.values()))))


# Write the intermediate tables with failed filters column but without removing anything
with pd.ExcelWriter(fname, mode='a') as writer:  
    df1new.to_excel(writer, sheet_name=sh1+'_all', na_rep = 'NA')
    df2new.to_excel(writer, sheet_name=sh2+'_all', na_rep = 'NA')


df1new_flt = df1new[(df1new.nsnp > 2) & (df1new.filter_samplesize == "TRUE") & (~df1new.outcome.isin(taxa_to_remove)) & (~df1new.exposure.isin(taxa_to_remove))]
df2new_flt = df2new[(df2new.nsnp > 2) & (df2new.filter_samplesize == "TRUE") & (~df2new.outcome.isin(taxa_to_remove)) & (~df2new.exposure.isin(taxa_to_remove))]

df1new_flt['BH_qval'] = multi.multipletests(df1new_flt['pval'], method = 'fdr_bh')[1]
df2new_flt['BH_qval'] = multi.multipletests(df2new_flt['pval'], method = 'fdr_bh')[1]

# Write the tables with basic filters and BH correction applied
with pd.ExcelWriter(fname, mode='a') as writer:  
    df1new_flt.to_excel(writer, sheet_name=sh1+'_flt', na_rep = 'NA')
    df2new_flt.to_excel(writer, sheet_name=sh2+'_flt', na_rep = 'NA')

