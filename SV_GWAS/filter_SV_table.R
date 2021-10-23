library(UpSetR)

run_qc_per_sv <- function(d, cohort_name, outpath, cr = 0.1){
    qc <- as.data.frame(t(apply(d, 2, function(x) table(factor(x, levels = c(0,1,NA)), exclude = NULL))))
    
    qc$num_called <- qc[,1] + qc[,2]
    qc$presence_rate <- qc[,2]/qc$num_called
    qc$call_rate <- qc$num_called / (qc[,1] + qc[,2] + qc[,3])

    pdf(paste0(outpath, ".qc_per_SV.CR",cr, ".pdf"))
    par(mfrow=c(2,2))
    hist(qc$num_called, breaks = 100, main = "Number of called dSV", col = "black")
    hist(qc$call_rate, breaks = 100, xlim = c(0,1),  main = "dSV call rate", col = "black")
    abline(v=cr, col = "red")
    plot(qc$presence_rate, qc$call_rate, pch = 16, cex = 0.5, col = "black", main = "Fraction of 1s vs call rate")
    abline(v=0.05, col = "red")
    abline(v=0.95, col = "red")
    abline(h=cr, col = "red")
    hist(qc$presence_rate, breaks = 100, xlim = c(0,1), col = "black", main = "dSV fraction of 1")
    abline(v=0.05, col = "red")
    abline(v=0.95, col = "red")
    dev.off()

    
    qc_05 <- qc[qc$num_called > 80 & qc$presence_rate > 0.05 & qc$presence_rate < 0.95 & qc$call_rate > cr,]
    #qc_05 <- qc[qc$call_rate > 0.15 & qc$presence_rate > 0.05 & qc$presence_rate < 0.95 & qc[,"1"] > 60 & qc[,"0"] > 60, ]

    cat(cohort_name, ": number SV before QC:", nrow(qc), "number SV after filtering:", nrow(qc_05), "\n", sep = " ")

    d <- d[colnames(d) %in% row.names(qc_05)]
    return(d)
}

run_qc_per_sample <- function(d, cohort_name, outpath){
    qc2 <-as.data.frame(apply(d, 1, function(y) length(which(! is.na(y)))))
    n2 <- ncol(d)
    colnames(qc2) <- "num_called"
    qc2$call_rate <- qc2$num_called/n2

    pdf(paste0(outpath, ".qc_per_sample.pdf"))
    hist(qc2$call_rate, breaks = 100, main = "Call rate per sample", col = "black")
    abline(v=0.05, col = "red")
    dev.off()
    qc2_05 <- qc2[qc2$call_rate > 0.05,]
    cat(cohort_name, ": number samples before QC:", nrow(qc2), "number samples after filtering:", nrow(qc2_05), "\n", sep = " ")

    d <- d[row.names(d) %in% row.names(qc2_05),]
    return(d)
}

setwd("/data/umcg-tifn/SV/SV_GWAS")
infile <- "/data/umcg-tifn/SV/profile/20210923_full_v3.0_10706samples/SV/20210923_full_deletionStructuralVariation_10706samples.tsv"
d <- read.delim(infile, header = T, sep = "\t", as.is = T, check.names = F, row.names = NULL)
cohorts <- read.delim("/data/umcg-tifn/SV/profile/20210923_full_v3.0_10706samples/SV/cleaned_file_list.tsv", header = T, sep = "\t", as.is = T, check.names = F)
cs <- c("LLD1", "300-OB", "500FG", "DAG3")
cohorts <- cohorts[cohorts$Cohort %in%  cs,]
sv_per_cohort <- matrix(nrow = ncol(d), ncol = length(cs))
row.names(sv_per_cohort) <- colnames(d)
colnames(sv_per_cohort) <- cs

for (c in cs){
    d_cohort <- d[d[,1] %in% cohorts[cohorts$Cohort == c, 1],]
    row.names(d_cohort) <- d_cohort[,1]
    d_cohort <- d_cohort[,2:ncol(d_cohort)]
    cr = 0.1
    if (c == "DAG3") cr = 0.05
    d_flt <- run_qc_per_sv(d_cohort, c, paste0("data/QC/", c, ".dSV_filtering"), cr)
    d_flt <- run_qc_per_sample(d_flt, c, paste0("data/QC/", c, ".dSV_filtering"))
    sv_per_cohort[row.names(sv_per_cohort) %in% colnames(d_flt), c] <- 1
}

#TODO 
sv_per_cohort <- as.data.frame(sv_per_cohort)
sv_per_cohort[is.na(sv_per_cohort)] <- 0
pdf("data/QC/dSV_overlap3.pdf")
upset(sv_per_cohort, order.by = "freq")
dev.off()

d[] <- lapply(as.data.frame(d), function(x) sub(1,2,x))
d[] <- lapply(as.data.frame(d), function(x) sub(0,1,x))
