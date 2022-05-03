args <- commandArgs(trailingOnly = TRUE)
library(UpSetR)
library(tibble)
library(dplyr)

infile <- args[1]
cohort_table <- args[2]

run_qc_per_sv <- function(d, cohort_name, outpath, cr = 0.1){
    qc <- as.data.frame(sapply(d, function(y) length(which(! is.na(y)))))
    
    n=nrow(d)
    nrow(qc)
    colnames(qc) <- "num_called"
    qc$call_rate <- qc[,1]/n
    
    pdf(paste0(outpath, ".qc_per_vSV.CR",cr, ".pdf"), useDingbats = F)
    par(mfrow=c(1,2))
    hist(qc$num_called, breaks = 100, main = "Number of called vSV", col = "black")
    hist(qc$call_rate, breaks = 100, xlim = c(0,1),  main = "vSV call rate", col = "black")
    abline(v=cr, col = "red")
    
    dev.off()

    qc_05 <- qc[qc$call_rate > cr,]
    
    cat(cohort_name, ": number SV before QC:", nrow(qc), "number SV after filtering:", nrow(qc_05), "\n", sep = " ")

    d <- d[colnames(d) %in% row.names(qc_05)]
    return(d)
}

run_qc_per_sample <- function(d, cohort_name, outpath){
    qc2 <-as.data.frame(apply(d, 1, function(y) length(which(! is.na(y)))))
    n2 <- ncol(d)
    colnames(qc2) <- "num_called"
    qc2$call_rate <- qc2[,1]/n2
    
    qc2_05 <- qc2[qc2$call_rate > 0.05,]
    
    pdf(paste0(outpath, "vSV.qc_per_sample.pdf"), useDingbats = F)
    hist(qc2$call_rate, breaks = 100, main = "Call rate per sample", col = "black")
    abline(v=0.05, col = "red")
    dev.off()
    qc2_05 <- qc2[qc2$call_rate > 0.05,]
    cat(cohort_name, ": number samples before QC:", nrow(qc2), "number samples after filtering:", nrow(qc2_05), "\n", sep = " ")

    d <- d[row.names(d) %in% row.names(qc2_05),]
    return(d)
}


setwd("/data/umcg-tifn/SV/SV_GWAS")
#infile <- "/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_full/20211212_full_deletionStructuralVariation_12388samples.tsv"
d <- read.delim(infile, header = T, sep = "\t", as.is = T, check.names = F, row.names = NULL)
conv <- read.delim("/data/umcg-tifn/SV/SV_GWAS/data/vsv_name_conversion_table.txt", header = T,  sep = "\t", as.is = T, check.names = F)
cohorts <- read.delim(cohort_table, header = T, sep = "\t", as.is = T, check.names = F)
cs <- c("LLD1", "300-OB", "500FG", "DAG3")
cohorts <- cohorts[cohorts$Cohort %in%  cs,]
sv_per_cohort <- matrix(nrow = ncol(d), ncol = length(cs))
row.names(sv_per_cohort) <- colnames(d)
colnames(sv_per_cohort) <- c("LLD", "300OB", "500FG", "DAG3")

for (c in cs){
    d_cohort <- d[d[,1] %in% cohorts[cohorts$Cohort == c, 1],]
    row.names(d_cohort) <- d_cohort[,1]
    d_cohort <- d_cohort[,2:ncol(d_cohort)]
    
    if (c == "LLD1") c = "LLD"
    if (c == "300-OB") c = "300OB"
    
    # get samples with genotypes
    geno <- read.delim(paste0("data/",c, ".covariates.txt") , header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    cat("Number of samples with vSVs: ", nrow(d_cohort), "\n")
    cat("Number of samples with genotypes and covariates: ", nrow(geno), "\n")
    sample_overlap <- row.names(d_cohort)[row.names(d_cohort) %in% row.names(geno)]
    cat ("Number overlapping samples: ", length(sample_overlap), "\n")
    d_cohort <- d_cohort[row.names(d_cohort) %in% sample_overlap,]
    
    #filter SVs
    d_flt <- run_qc_per_sv(d_cohort, c, paste0("data/QC/", c, ".vSV_filtering"))
    d_flt <- run_qc_per_sample(d_flt, c, paste0("data/QC/", c, ".vSV_filtering"))

    sv_per_cohort[row.names(sv_per_cohort) %in% colnames(d_flt), c] <- 1

    # normalize: apply INT
    d_norm <- apply(d_flt, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
    

    # rename SVs
    new_sv_ids <-(conv[match(colnames(d_norm), conv$sv_id, nomatch = 0),"new_sv_id"])
    length(new_sv_ids) == ncol(d_norm)
    colnames(d_norm) <- new_sv_ids

    d_norm <- rownames_to_column(as.data.frame(d_norm),var = "#IID")
    #write filtered table
    write.table(d_norm, file = paste0("data/", c, ".vSV.filtered.txt"), sep = "\t", quote = F, row.names = F) 
}


# plot number of SVs and their overlap between cohorts
sv_per_cohort <- as.data.frame(sv_per_cohort)
sv_per_cohort[is.na(sv_per_cohort)] <- 0
pdf("data/QC/vSV_overlap.pdf", useDingbats = F)
upset(sv_per_cohort, order.by = "freq")
dev.off()


sv_per_cohort2 <- as.data.frame(sv_per_cohort[rowSums(sv_per_cohort, na.rm = T) > 1,])
new_sv_ids2 <-(conv[match(row.names(sv_per_cohort2), conv$sv_id, nomatch = 0),"new_sv_id"])
length(new_sv_ids2) == nrow(sv_per_cohort2)
row.names(sv_per_cohort2) <- new_sv_ids2
sv_per_cohort2$cohorts <- NA

sv_per_cohort2$cohorts <- apply(sv_per_cohort2, 1, function(x) {paste(colnames(sv_per_cohort2)[which(x==1)], collapse = ",")})

write.table(sv_per_cohort2, file = "data/vSV_per_cohort.txt", sep = "\t", quote = F, col.names = NA) 


# SV in any number of cohorts
sv_per_cohort <- sv_per_cohort[2:nrow(sv_per_cohort),]
new_sv_ids3 <-(conv[match(row.names(sv_per_cohort), conv$sv_id, nomatch = 0),"new_sv_id"])
length(new_sv_ids3) == nrow(sv_per_cohort)
row.names(sv_per_cohort) <- new_sv_ids3
sv_per_cohort <- as.data.frame(sv_per_cohort)
sv_per_cohort$cohorts <- NA
sv_per_cohort$cohorts <- apply(sv_per_cohort, 1, function(x) {paste(colnames(sv_per_cohort)[which(x==1)], collapse = ",")})
sv_per_cohort <- sv_per_cohort[sv_per_cohort$cohorts != "",]

write.table(sv_per_cohort, file = "data/vSV_per_cohort_all_cohorts.txt", sep = "\t", quote = F, col.names = NA) 
