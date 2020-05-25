args <- commandArgs(trailingOnly = TRUE)

library(TwoSampleMR)
library(MRInstruments)

source("/groups/umcg-lld/tmp03/umcg-dzhernakova/umcg_scripts/MR/run_MR.R")
access_token_fname = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/v2/mrbase.oauth"
out_path = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/BCAA/300OB/MR_300OB_pQTLs.FDR0.05-AA.results.txt"
qtl_fname = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/BCAA/300OB/300ob_eQTLsFDR0.05-ProbeLevel.forMR.txt"
outcome_ids = c("897", "873", "940")


pqtl_table <- read.table(qtl_fname, header = T, sep = "\t", as.is = T, check.names = F)
res_table <- data.frame()

for (outcome_id in outcome_ids){
    for (p in unique(pqtl_table$Phenotype){
        paste(outcome_id, p) 
        out_dat = NULL
        exp_dat = NULL
        qtls_subs <- pqtl_table[pqtl_table$Phenotype == p,]

        tryCatch({
            out_dat <- extract_outcome_data(snps = qtls_subs$SNP, outcomes = outcome_id, access_token = access_token_fname)
            exp_dat <- clump_data(format_data(qtls_subs, snps = out_dat$SNP, type = "exposure"))
        }, error=function(e) NULL)

        if (!is.null(out_dat) & !is.null(exp_dat)){
            res <- run_mr(exp_dat, out_dat, exp_table = qtls_subs)
            if (!is.null(res)){
                res_table <- rbind(res_table, res)
            }
        }
    }
}
write.table(res_table, file =  out_path , sep = "\t", quote = F, col.names = NA)

print("Finished")
