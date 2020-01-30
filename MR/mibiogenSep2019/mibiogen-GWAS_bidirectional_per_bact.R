args <- commandArgs(trailingOnly = TRUE)

library(TwoSampleMR)
library(MRInstruments)

source("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/v2/run_MR.R")

ao <- available_outcomes(access_token = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/v2/mrbase.oauth")

outcome_ids <- c("2","7","9","10","12","22","24","30","31","32","44","45","54","72","276","292","294","297","298","798","800","801","802","803","805","806","812","814","815","833","970","990","996","1024","1025","1029","1054","1058","1085","1086","1090","1108","1109","1110","1111","1112","302", "299", "300", "301")

mibiogen_path <- args[1]
bact <- basename(mibiogen_path)

pval_thres = 5e-08

out_filebase <- paste0("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/mibiogenSep2019/per_bact3/", pval_thres)

print(paste0("Taking mbQTLs from: ", mibiogen_path))
print(paste0("p-value threshold for mbQTLs when used as exposure: ",pval_thres)) 
print(paste0("Writing to this folder:",out_filebase))
phenotype="gwas"
print(phenotype)

res_table <- data.frame()
mibiogen_table <- read.table(gzfile(mibiogen_path), header = T, sep = "\t", as.is = T, check.names = F)
#
# GWAS - mibiogen
#
print ("GWAS -> mibiogen")
#for (exp_id in outcome_ids){
#  print(paste0(exp_id, " -> ", bact))
#  out_dat = NULL
#  exp_dat = NULL
#  exp_dat0 = NULL
#  tryCatch({
#    exp_dat0 <- extract_instruments(outcomes=exp_id, clump = FALSE, access_token = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/v2/mrbase.oauth")
#    if (nrow(exp_dat0) > 0){
#     out_dat <- format_data(mibiogen_table, snps = exp_dat0$SNP, type = "outcome")
#     exp_dat <- clump_data(exp_dat0[exp_dat0$SNP %in% out_dat$SNP,])
#    }
#  },error=function(e) {
#    message(paste("NB!! Data extraction didn't work for", exp_id, " -> ", bact))
#  })  
#  if (!is.null(out_dat) & !is.null(exp_dat)){
#    res <- run_mr(exp_dat, out_dat)
#    if (!is.null(res)){
#      res_table <- rbind(res_table, res)
#    }
#  }
#}
#
#write.table(res_table, file =  paste0(out_filebase, "/", phenotype, "-mibiogen.", pval_thres, ".", bact) , sep = "\t", quote = F, col.names = NA)
#


#
# mibiogen -> GWAS with a less stringent threshold 1e-5
#
pval_thres = 5e-8
out_filebase <- paste0("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/mibiogenSep2019/per_bact3/", pval_thres)

print(paste0("Taking mbQTLs from: ", mibiogen_path))
print(paste0("p-value threshold for mbQTLs when used as exposure: ",pval_thres))
res_table <- data.frame()

exp_table <- mibiogen_table[mibiogen_table$pval < as.numeric(pval_thres),]
  if (nrow(exp_table) > 0){
    for (outcome_id in outcome_ids){
      print(paste0(bact, " -> ", outcome_id))
      out_dat = NULL
      exp_dat = NULL
      tryCatch({
        out_dat <- extract_outcome_data(snps = exp_table$SNP, outcomes = outcome_id, access_token = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/v2/mrbase.oauth")
        exp_dat <- clump_data(format_data(exp_table, snps = out_dat$SNP, type = "exposure"))
      }, error=function(e) {
        message(paste("NB!! Data extraction didn't work for", bact, " -> ", outcome_id))
      }) 
      if (!is.null(out_dat) & !is.null(exp_dat)){
        res <- run_mr(exp_dat, out_dat)
        if (!is.null(res)){
          res_table <- rbind(res_table, res)
        }
      }
    }
  }

write.table(res_table, file = paste0(out_filebase, "/mibiogen-", phenotype, ".", pval_thres, ".", bact), sep = "\t", quote = F, col.names = NA)