f = open("C:/Users/Dasha/work/UMCG/data/MR/results2/mibiogen/mibiogenSep2019/tmp_res4filtering.txt")
out = open("C:/Users/Dasha/work/UMCG/data/MR/results2/mibiogen/mibiogenSep2019/tmp_res4filtering_res.txt", 'w')
header = f.readline().rstrip()
cols = {c : i for i, c in enumerate(header.split('\t'))}
print cols
out.write(header + "\tfailed_filters" + "\n")
for l in f:
    spl = l.rstrip().split("\t")
    flt_line = '-'
    if int(spl[cols['nsnp']]) < 3:
        flt_line += ';Number of SNPs'
    if not spl[cols['egger_intercept_pval']] == 'NA':
        if float(spl[cols['egger_intercept_pval']]) < 0.05:
            flt_line += ';Egger intercept p-value'
    if not spl[cols['weighted_median_pval']] == 'NA':        
        if float(spl[cols['weighted_median_pval']]) > 0.05:
            flt_line += ';Weighted median p-value'
    if not spl[cols['mr_presso_global']] == 'NA':
        if not spl[cols['mr_presso_global']] == '<1e-04':
            if not spl[cols['mr_presso_outlier_cor_pval']] == 'NA':
                if float(spl[cols['mr_presso_global']]) < 0.05 and float(spl[cols['mr_presso_outlier_cor_pval']]) > 0.05:
                    flt_line += ';MR PRESSO outlier test'
    if not spl[cols['mr_presso_pval']] == 'NA':
        if float(spl[cols['mr_presso_pval']]) > 0.05:
            flt_line += ';MR PRESSO p-value'
    if spl[cols['filter_loo']] == 'FALSE':
        flt_line += ';Leave-one-out analysis'
    out.write(l.rstrip() + "\t" + flt_line.replace("-;", "", 1) + "\n")
out.close()