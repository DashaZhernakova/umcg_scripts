
### Age and gender
pheno <- read.delim("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v1/tab_separated_labels/Participant (Pa_99_G).dat", header = T, sep = "\t", as.is = T, check.names = F)
row.names(pheno) <- pheno$PSEUDOIDEXT
pheno <- subset(pheno, select = c(GESLACHT, AGE_1A1, AGE_2A1))
colnames(pheno) <- c("gender_F1M2", "age", "age_2a")
pheno$gender_F1M2 <- sub("Female", 1, pheno$gender_F1M2)
pheno$gender_F1M2 <- sub("Male", 2, pheno$gender_F1M2)


format_phenotypes <- function(dtable, phenos2select){
    dtable_a1 <- dtable[dtable$ENCOUNTERCODE == "Baseline assessment (1A)",]
    dtable_a2 <- dtable[dtable$ENCOUNTERCODE == "Second assessment (2A)",]
    row.names(dtable_a1) <- dtable_a1$PSEUDOIDEXT
    row.names(dtable_a2) <- dtable_a2$PSEUDOIDEXT
    dtable_a1 <- as.data.frame(dtable_a1[,phenos2select, drop = F])
    dtable_a2 <- as.data.frame(dtable_a2[,phenos2select, drop = F])
    return (list(dtable_a1, dtable_a2))
}

### Measurements
blood <- read.delim("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v1/tab_separated_labels/Laboratory_assessment_Blood.dat", header = T, sep = "\t", as.is = T, check.names = F)
blood <- blood[blood$NUCHTER == "Yes", ]
blood2select <- c("GLU", "HB1C", "CHO", "HDC", "LDC", "TGL")
new_blood <- format_phenotypes(blood, blood2select)
blood_a1 <- new_blood[[1]]
blood_a2 <- new_blood[[2]]


bp <- read.delim("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v1/tab_separated_labels/Measurement_Blood_pressure.dat", header = T, sep = "\t", as.is = T, check.names = F)
bp2select <- c("SBP")
new_bp <- format_phenotypes(bp, bp2select)
bp_a1 <- new_bp[[1]]
bp_a2 <- new_bp[[2]]

anthro <- read.delim("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v1/tab_separated_labels/Measurement_Anthropometry.dat", header = T, sep = "\t", as.is = T, check.names = F)
anthro2select <- c("BMI")
new_anthro <- format_phenotypes(anthro, anthro2select)
anthro_a1 <- new_anthro[[1]]
anthro_a2 <- new_anthro[[2]]

# add smoking

### merge phenotypes
ida1 <- intersect(intersect(row.names(pheno), row.names(blood_a1)), intersect(row.names(bp_a1), row.names(anthro_a1)))
ida2 <- intersect(intersect(row.names(pheno), row.names(blood_a2)), intersect(row.names(bp_a2), row.names(anthro_a2)))
ids <- intersect(ida1, ida2)
length(ida1)
length(ida2)
length(ids)

a1 <- cbind(pheno[ida1,], blood_a1[ida1,], bp_a1[ida1,], anthro_a1[ida1,])
colnames(a1) <- c(colnames(pheno), colnames(blood_a1), colnames(bp_a1), colnames(anthro_a1))

a2 <- cbind(pheno[ida2,], blood_a2[ida2,], bp_a2[ida2,], anthro_a2[ida2,])
colnames(a2) <- c(colnames(pheno), colnames(blood_a2), colnames(bp_a2), colnames(anthro_a2))

m <- cbind(pheno[ids,], blood_a1[ids,], bp_a1[ids,], anthro_a1[ids,], blood_a2[ids,], bp_a2[ids,], anthro_a2[ids,])
colnames(m) <- c(colnames(pheno), paste0(c(colnames(blood_a1), colnames(bp_a1), colnames(anthro_a1)), "_1a"), paste0(c(colnames(blood_a2), colnames(bp_a2), colnames(anthro_a2)), "_2a"))


m$age_bin <- cut(m$age, breaks = c(0, seq(20, 80, 5), 100))
phenos <- c(blood2select, anthro2select, bp2select)
res_table <- data.frame(matrix(ncol = 7, nrow = length(unique(m$gender_F1M2)) * length(unique(m$age_bin)) * length(phenos) ))
colnames(res_table) <- c("pheno", "gender_F1M2", "age_bin", "n_samples", "pvalue", "mean_dif", "tvalue")
i <- 1
for (gender in unique(m$gender_F1M2)) {
    for (pheno in phenos){
        for (a in unique(m$age_bin)) {
            x <- m[m$age_bin == a & m$gender_F1M2 == gender, c(paste0(pheno,"_1a"), paste0(pheno,"_2a"))]
            x <- na.omit(x)
            res <- t.test(x[,2], x[,1], paired = TRUE)
            #cat (pheno, a, nrow(x), res$p.value, res$estimate, res$statistic, "\n", sep = ", ")
            res_table[i,] <- c(pheno, gender, a, nrow(x), res$p.value, res$estimate, res$statistic)
            i <- i + 1
        }
    }
}

res_table <- res_table[res_table$age_bin != "(0,20]" & res_table$age_bin != "(80,100]", ]

##
### Diseases
##

# now from LL v2
idconv <- read.delim("/groups/umcg-lifelines/prm03/releases/pheno_lifelines/v1/phenotype_linkage_file_project_pseudo_id.txt", header = T, sep = "\t", as.is = T, check.names = F)

d1a <- read.delim("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v2/1a_q_1_cvd.txt", header = T, sep = "\t", as.is = T, check.names = F)
d1a <- d1a[,c("PROJECT_PSEUDO_ID","aneurysm_diagnosis_adu_q_1","angioplasty_bypass_adu_q_1","arrhythmia_presence_adu_q_1",  "carotid_stenosis_adu_q_1","heartattack_presence_adu_q_1","heartfailure_presence_adu_q_1",  "highcholesterol_presence_adu_q_1","stroke_presence_adu_q_1")]
healthy = d1a[which(apply(d1a[,seq(4,ncol(d1a))],1,function(x) all(x == 0))),"PROJECT_PSEUDO_ID"]

d1b <- read.delim("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v2/1b_q_1_cvd.txt", header = T, sep = "\t", as.is = T, check.names = F)
new_cvd_1b <- d1b[which(apply(d1b[,seq(2,ncol(d1b))],1,function(x) any(x == 1))),"PROJECT_PSEUDO_ID"]

d1c <- read.delim("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v2/1c_q_1_cvd.txt", header = T, sep = "\t", as.is = T, check.names = F)
new_cvd_1c <- d1c[which(apply(d1c[,seq(2,ncol(d1c))],1,function(x) any(x == 1))),"PROJECT_PSEUDO_ID"]

d2 <- read.delim("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v2/2a_q_1_cvd.txt", header = T, sep = "\t", as.is = T, check.names = F)
new_cvd_2a <- d2[which(apply(d2[,seq(2,ncol(d2))],1,function(x) any(x == 1))),"PROJECT_PSEUDO_ID"]

d3a <- read.delim("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v2/2a_q_1_cvd.txt", header = T, sep = "\t", as.is = T, check.names = F)
new_cvd_3a <- d3a[which(apply(d3a[,seq(2,ncol(d3a))],1,function(x) any(x == 1))),"PROJECT_PSEUDO_ID"]

