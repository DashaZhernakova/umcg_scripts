setwd("/groups/umcg-lld/scr01/dasha/MR/results/AA_T2D/")
source("/groups/umcg-lld/scr01/dasha/MR/results/run_MR.R")

out_filebase <- "MR_AA-pQTL.tested_subset.1e-7.results"
#qtl_fname <- "/groups/umcg-lld/scr01/dasha/MR/results/AA_T2D/all_cis+trans-pQTLs.for_MR.rs.txt"
qtl_fname <- "/groups/umcg-lld/scr01/dasha/MR/results/AA_T2D/pQTLs.BCAA_SNPs.txt"

pqtl_table <- read.table(qtl_fname, header = T, sep = "\t", as.is = T, check.names = F)
prots <- c("CPA1","JAM-A","uPA","CPB1","CHI3L1","t-PA","SCGB3A2","EGFR","CD93","IL-18BP","PON3","CTSZ","MMP-3","RARRES2","ICAM-2","KLK6","IGFBP-2","PECAM-1","NT-pro-BNP","CCL16")

aa_ids <- c("897", "873", "940")

res_table_forw <- data.frame()
res_table_rev <- data.frame()

for (aa_id in aa_ids){
  for (p in prots){
    # AA -> pQTLs
    paste(aa_id, p) 
    out_dat = NULL
    exp_dat = NULL
    qtls_subs <- pqtl_table[pqtl_table$Phenotype == p,]
    tryCatch({
      exp_dat <- extract_instruments(outcomes=aa_id, p1=1e-7, p2=1e-7,access_token = "/groups/umcg-lld/scr01/dasha/MR/results/mrbase.oauth")
      out_dat <- format_data(qtls_subs, snps = exp_dat$SNP, type = "outcome")
    }, error=function(e) NULL)
    if (!is.null(out_dat) & !is.null(exp_dat)){
      res <- run_mr(exp_dat, out_dat)
      if (!is.null(res)){
        res_table_forw <- rbind(res_table_forw, res)
      }
    }
    
  }
}

res_table_forw$BH_qval = NA
flt <- which((res_table_forw$nsnp > 2 | "cis" %in% res_table_forw$type) & res_table_forw$egger_intercept_pval > 0.05 & res_table_forw$heterogeneity_Q_pval > 0.05)
res_table_forw[flt, "BH_qval"] <- p.adjust(res_table_forw[flt, "pval"], method = "BH")
write.table(res_table_forw, file = paste0(out_filebase, ".txt") , sep = "\t", quote = F, col.names = NA)
