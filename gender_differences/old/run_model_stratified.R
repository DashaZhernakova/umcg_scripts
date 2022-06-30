args <- commandArgs(trailingOnly = TRUE)
source("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/scripts/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
setwd("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/factors")

rm_outliers <- function(merged_tab){
    merged_tab <- na.omit(merged_tab)
    merged_tab <- merged_tab[(merged_tab$age < 80) & (merged_tab$age >= 20),]
  
  w <- merged_tab[merged_tab$gender_F1M2 == 1,]
  m <- merged_tab[merged_tab$gender_F1M2 == 2,]
  mq1 <- quantile(m[,1], probs = 0.25)
  mq3 <- quantile(m[,1], probs = 0.75)
  miqr <- mq3 - mq1
  m_clean <- m[m[,1] < mq3 + 1.5*miqr & m[,1] > mq1 - 1.5*miqr,]
  
  wq1 <- quantile(w[,1], probs = 0.25)
  wq3 <- quantile(w[,1], probs = 0.75)
  wiqr <- wq3 - wq1
  w_clean <- w[w[,1] < wq3 + 1.5*wiqr & w[,1] > wq1 - 1.5*wiqr,]
  
  tab_nooutliers <- rbind(w_clean, m_clean)
  return(tab_nooutliers)
}



traits_path <- "../v4/data/LL_phenotypes_merged_all.log_some.v5.txt"
pheno_path <- "factors_final.txt"
#pheno="CHO"
#med="statins"

phenos <- c("AF", "ALT", "AST","CA", "FOS", "NAA", "GLU",  "DBP", "CHO", "HDC", "LDC", "TGL", "HALB", "GR", "ER",  "HT", "BALB")
meds <- c("AF_med", "ALT_med", "AST_med","calcium", "FOS_med", "NAA_med", "T2D_med",  "antihypertensives", "statins", "statins", "statins", "statins", "statins", "GR_med", "ER_med", "HT_med", NA)


out_path<-"/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/factors/results/"
# read phenotype traits of interest
traits <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

# read age, gender and other covariate phenotypes
pheno <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)


#order samples in the two tables
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ), , drop = F]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0), ]
all(row.names(traits_m) == row.names(pheno_m))
num_traits <- ncol(traits_m)

cnt = 1
pheno_num = 1
for (pheno in phenos){
    med <- meds[pheno_num]

    if (!is.na(med)){
        merged_tab <- cbind(traits_m[, pheno], pheno_m[,1:9], pheno_m[,med])
        colnames(merged_tab)[1] <- "phenotype"
        colnames(merged_tab)[ncol(merged_tab)] <- "med"

        merged_tab <- rm_outliers(merged_tab)

        merged_tab <- mutate(merged_tab, med = ordered(med, levels = c('0', '1')))

        base_formula <- as.formula("phenotype ~ smoking + med + s(stress_1y, k = 7) + s(stress_chronic, k = 7) +
        s(phys_activity_total) + s(phys_activity_intensive) + s(diet) +
        s(alcohol)")
    } else {
        merged_tab <- cbind(traits_m[, pheno], pheno_m[,1:9])
        colnames(merged_tab)[1] <- "phenotype"
        
        merged_tab <- rm_outliers(merged_tab)
        
        base_formula <- as.formula("phenotype ~ smoking  + s(stress_1y, k = 7) + s(stress_chronic, k = 7) +
        s(phys_activity_total) + s(phys_activity_intensive) + s(diet) +
        s(alcohol)")
    }
    cat("\n\n-----------", pheno, "-----------\n\n")
    
    merged_tab <- mutate(merged_tab, gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
    merged_tab$smoking <- gsub("2", "0", merged_tab$smoking)
    merged_tab <- mutate(merged_tab, smoking = ordered(smoking, levels = c('0', '1')))

    samplesizes <- c()
    for (min_age in seq(20,65,15)){
      #cat(min_age, sex, "\n")
      max_age = min_age + 15
      w_n <- nrow(merged_tab[merged_tab$gender_F1M2 == "1" & merged_tab$age > min_age & merged_tab$age < max_age,])
      m_n <- nrow(merged_tab[merged_tab$gender_F1M2 == "2" & merged_tab$age > min_age & merged_tab$age < max_age,])
      samplesizes <- c(samplesizes, w_n, m_n)
      cat(min_age, w_n, m_n, "\n", sep = " ")
    }
    nsamples <- min(samplesizes)
    
    for (min_age in seq(20,65,15)){
        #cat(min_age, sex, "\n")
        max_age = min_age + 15
        w <- merged_tab[merged_tab$gender_F1M2 == "1" & merged_tab$age > min_age & merged_tab$age < max_age,]
        m <- merged_tab[merged_tab$gender_F1M2 == "2" & merged_tab$age > min_age & merged_tab$age < max_age,]
        
        # make equal samplesizes for p-value comparison 
        # n_w <- nrow(w)
        # n_m <- nrow(m)
        # if (n_w < n_m){
        #     m <- m[sample(n_m, n_w),]
        # } else if (n_m < n_w){
        #     w <- w[sample(n_w, n_m),]
        # }
        w <- w[sample(nrow(w), nsamples), ]
        m <- m[sample(nrow(m), nsamples), ]
        
        # women
        gam_fit <- gam(base_formula, data = w, method = "REML", select=T)
        s <- summary(gam_fit)
        p_table <- s$p.table
        s_table <- s$s.table

        rows <- c(row.names(p_table), row.names(s_table))
        res_line <- c(p_table[,4], s_table[,4])
        if (cnt == 1){
            res <- data.frame(matrix(nrow = length(rows), ncol = length(phenos)*8))
            row.names(res) <- rows
        }
        colnames(res)[cnt] = paste(pheno, min_age, "w", sep = "_")
        #res_m <- match(row.names(res), rows, nomatch = 0)
        res[match(rows,row.names(res), nomatch = 0),cnt] <- res_line[match(row.names(res), rows, nomatch = 0)]
        cnt <- cnt + 1

        # men
        gam_fit <- gam(base_formula, data = m, method = "REML", select=T)
        s <- summary(gam_fit)
        p_table <- s$p.table
        s_table <- s$s.table

        rows <- c(row.names(p_table), row.names(s_table))
        res_line <- c(p_table[,4], s_table[,4])
        colnames(res)[cnt] = paste(pheno, min_age, "m", sep = "_")
        #m <- match(row.names(res), rows, nomatch = 0)
        res[match(rows,row.names(res), nomatch = 0),cnt] <- res_line[match(row.names(res), rows, nomatch = 0)]
        cnt <- cnt + 1
    }
    pheno_num = pheno_num + 1
}
write.table(res, file = paste0(out_path, "/stratified_gam_main_sameN.txt"), sep = "\t", quote = F, col.names = NA)



library(corrplot)
pdf(paste0(out_path, "/lifestyle_factors_stratified_corrplot2_sameN2.pdf"), width = 20, height = 5, useDingbats=FALSE)
cuts <- apply(res, 2, function(x) {3-as.numeric(cut(x, c(-Inf,0.001, 0.05,1), labels=0:2))})

row.names(cuts) <- row.names(res)
corrplot(cuts,is.corr = F, tl.col = "black", tl.cex = 0.7)
dev.off()
