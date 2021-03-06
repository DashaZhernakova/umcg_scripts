args <- commandArgs(trailingOnly = TRUE)

library(TwoSampleMR)
library(MRInstruments)
library(MRPRESSO)

source("/groups/umcg-lld/tmp03/umcg-dzhernakova/umcg_scripts/MR/run_MR.R")


mibiogen_path <- args[1]
bact <- basename(mibiogen_path)

pval_thres = 5e-08
#pval_thres =  1e-05


out_filebase <- paste0("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/mibiogenOct2020/per_bact/", pval_thres)

print(paste0("Taking mbQTLs from: ", mibiogen_path))
print(paste0("p-value threshold for mbQTLs when used as exposure: ",pval_thres))
print(paste0("Writing to this folder:",out_filebase))

res_table <- data.frame()

file.names <- dir("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/data/UKB/original_files_Neale/UKB_formatted/", pattern =".for_MR.txt.gz")

mibiogen_table <- read.table(gzfile(mibiogen_path), header = T, sep = "\t", as.is = T, check.names = F)
exp_table <- mibiogen_table[mibiogen_table$pval < as.numeric(pval_thres),]
#for(i in 1:length(file.names)){
#  print(file.names[i])
#  pheno_table <- read.table(gzfile(paste0("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/data/UKB/",file.names[i])), header = T, sep = "\t", as.is = T, check.names = F)
#  #
#  # mibiogen -> GWAS
#  #
#  if (nrow(exp_table) > 0){
#      out_dat = NULL
#      exp_dat = NULL
#      tryCatch({
#        out_dat <- format_data(pheno_table, snps = exp_table$SNP, type = "outcome")
#        exp_dat <- clump_data(format_data(exp_table, snps = out_dat$SNP, type = "exposure"))
#      },error=function(e) NULL)
#      if (!is.null(out_dat) & !is.null(exp_dat)){
#        res <- run_mr(exp_dat, out_dat)
#        if (!is.null(res)){
#          res_table <- rbind(res_table, res)
#        }
#      }
#  }
#}
#res_table$filter_before_BH <- (res_table$nsnp > 2 | "cis" %in% res_table$type) & res_table$egger_intercept_pval > 0.05 & res_table$heterogeneity_Q_pval > 0.05
#write.table(res_table, file = paste0(out_filebase, "/mibiogen-UKB.Fstat.", pval_thres, ".", bact), sep = "\t", quote = F, col.names = NA)

res_table <- data.frame()

  #
  # GWAS - mibiogen
  #
  print ("GWAS -> mibiogen")
for(i in 1:length(file.names)){
  print(file.names[i])
    pheno_table <- read.table(gzfile(paste0("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/data/UKB/original_files_Neale/UKB_formatted/",file.names[i])), header = T, sep = "\t", as.is = T, check.names = F)
    pheno_table_subs <- pheno_table[pheno_table$pval < 5e-08,]
    out_dat = NULL
    exp_dat = NULL
    tryCatch({
      out_dat <- format_data(mibiogen_table, snps = pheno_table_subs$SNP, type = "outcome")
      exp_dat <- clump_data(format_data(pheno_table_subs, snps = out_dat$SNP, type = "exposure"))
    },error=function(e) NULL)
    if (!is.null(out_dat) & !is.null(exp_dat)){
      res <- run_mr(exp_dat, out_dat)
      if (!is.null(res)){
        res_table <- rbind(res_table, res)
      }
    }
  }


#res_table$BH_qval <- p.adjust(res_table$pval, method = "BH")
#res_table$filter_before_BH <- (res_table$nsnp > 2 | "cis" %in% res_table$type) & res_table$egger_intercept_pval > 0.05 & res_table$heterogeneity_Q_pval > 0.05

write.table(res_table, file =  paste0(out_filebase, "/UKB-mibiogen.", pval_thres, ".", bact) , sep = "\t", quote = F, col.names = NA)


#
# mibiogen -> GWAS with a less stringent threshold 1e-5
#
pval_thres = 1e-5
out_filebase <- paste0("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/mibiogenOct2020/per_bact/", pval_thres)

print(paste0("Taking mbQTLs from: ", mibiogen_path))
print(paste0("p-value threshold for mbQTLs when used as exposure: ",pval_thres))
res_table <- data.frame()

exp_table <- mibiogen_table[mibiogen_table$pval < as.numeric(pval_thres),]

for(i in 1:length(file.names)){
  print(file.names[i])
  pheno_table <- read.table(gzfile(paste0("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/data/UKB/original_files_Neale/UKB_formatted/",file.names[i])), header = T, sep = "\t", as.is = T, check.names = F)
  #
  # mibiogen -> GWAS
  #
  if (nrow(exp_table) > 0){
      out_dat = NULL
      exp_dat = NULL
      tryCatch({
        out_dat <- format_data(pheno_table, snps = exp_table$SNP, type = "outcome")
        exp_dat <- clump_data(format_data(exp_table, snps = out_dat$SNP, type = "exposure"))
      },error=function(e) NULL)
      if (!is.null(out_dat) & !is.null(exp_dat)){
        res <- run_mr(exp_dat, out_dat)
        if (!is.null(res)){
          res_table <- rbind(res_table, res)
        }
      }
  }
}

#res_table$filter_before_BH <- (res_table$nsnp > 2 | "cis" %in% res_table$type) & res_table$egger_intercept_pval > 0.05 & res_table$heterogeneity_Q_pval > 0.05
write.table(res_table, file = paste0(out_filebase, "/mibiogen-UKB.", pval_thres, ".", bact), sep = "\t", quote = F, col.names = NA)
