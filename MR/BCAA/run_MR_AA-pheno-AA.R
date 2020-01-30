setwd("/groups/umcg-lld/scr01/dasha/MR/results/AA_T2D/")
source("/groups/umcg-lld/scr01/dasha/MR/results/run_MR.R")

out_filebase <- "MR_AA-lipids_results"

ao <- available_outcomes(access_token = "/groups/umcg-lld/scr02/dasha/MR/results/v2/mrbase.oauth")
#lipids
t2d_ids <- c(ao[ao$subcategory == "Lipid" & ao$author == "Kettunen", "id"])
#intrinsic factors
#Creatinine, Glucose, HbA1c, Total cholesterol, HDL-C, LDL-C, TG, Citrullin, Haemoglobin, Insulin, Lymph%
#Leptin, adiponectin
#t2d_ids <- c("309", "850", "UKB-a:333", "1099", "1103", "1104", "1105", "18", 
#	"418", "422", "756", "757", "772", "773", "776", "777", "859",
#	"758",
#	"933", "301", "782",
#	"299", "780",
#	"300", "781",
#	"934", "302", "783",
#	"356",
#	"270", "271", "272",
#	"761", "762", "768", "774", "775", "778", "779", "767",
#	"119", "125", "134", "139", "140",
#	"1002", "1003",
#	"1")


# Amino acids from Kettunen et al
#aa_ids <- c("840","850","860","866","873","897","919","938","939","940")
# all BCAA 
aa_ids <- ao[ao$trait %in% c("Leucine", "Isoleucine", "Valine"), "id"]


res_table_forw <- data.frame()
res_table_rev <- data.frame()

for (aa_id in aa_ids){
  for (t2d_id in t2d_ids){
    # AA -> T2D
    tryCatch({
    exp_dat <- extract_instruments(outcomes = aa_id, access_token = "/groups/umcg-lld/scr02/dasha/MR/results/v2/mrbase.oauth")
    out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = t2d_id, access_token = "/groups/umcg-lld/scr02/dasha/MR/results/v2/mrbase.oauth")
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
flt <- which(res_table_forw$nsnp > 2 & res_table_forw$egger_intercept_pval > 0.05 & res_table_forw$heterogeneity_Q_pval > 0.05)
res_table_forw[flt, "BH_qval"] <- p.adjust(res_table_forw[flt, "pval"], method = "BH")
write.table(res_table_forw, file = paste0(out_filebase, ".forward.txt") , sep = "\t", quote = F, col.names = NA)



for (aa_id in aa_ids){
  for (t2d_id in t2d_ids){
    # T2D -> AA
    tryCatch({
    exp_dat <- extract_instruments(outcomes = t2d_id, access_token = "/groups/umcg-lld/scr02/dasha/MR/results/v2/mrbase.oauth")
    out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = aa_id, access_token = "/groups/umcg-lld/scr02/dasha/MR/results/v2/mrbase.oauth")
    }, error=function(e) NULL) 
    if (!is.null(out_dat) & !is.null(exp_dat)){
      res <- run_mr(exp_dat, out_dat)
      if (!is.null(res)){
        res_table_rev <- rbind(res_table_rev, res)
      }
    }
  }
}



res_table_rev$BH_qval = NA
flt <- which(res_table_rev$nsnp > 2 & res_table_rev$egger_intercept_pval > 0.05 & res_table_rev$heterogeneity_Q_pval > 0.05)
res_table_rev[flt, "BH_qval"] <- p.adjust(res_table_rev[flt, "pval"], method = "BH")
write.table(res_table_rev, file = paste0(out_filebase, ".reverse.txt") , sep = "\t", quote = F, col.names = NA)