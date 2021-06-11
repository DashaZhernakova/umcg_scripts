d <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/1a_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")
d[d == "$6"] <- 0
d[d == "3"] <- NA
d[d == ""] <- NA
d[d == "2"] <- 0

d$arrhythmia_presence_adu_q_1 <- ifelse(is.na(d$arrhythmia_presence_adu_q_1), d$arrhythmia_presence_adu_q_2, d$arrhythmia_presence_adu_q_1)
d$heartattack_presence_adu_q_1 <- ifelse(is.na(d$heartattack_presence_adu_q_1), d$heartattack_presence_adu_q_2, d$heartattack_presence_adu_q_1)
d$heartfailure_presence_adu_q_1 <- ifelse(is.na(d$heartfailure_presence_adu_q_1), d$heartfailure_presence_adu_q_2, d$heartfailure_presence_adu_q_1)
d$stroke_presence_adu_q_1 <- ifelse(is.na(d$stroke_presence_adu_q_1), d$stroke_presence_adu_q_2, d$stroke_presence_adu_q_1)

d1a <- d[,c("PROJECT_PSEUDO_ID","AGE","GENDER","aneurysm_diagnosis_adu_q_1","aneurysm_presence_adu_q_2","angioplasty_bypass_adu_q_1","arrhythmia_presence_adu_q_1", "atherosclerosis_presence_adu_q_1","carotid_stenosis_adu_q_1","heartattack_presence_adu_q_1","heartfailure_presence_adu_q_1","heartvalve_presence_adu_q_1","highcholesterol_presence_adu_q_1","stenosis_presence_adu_q_1","stroke_presence_adu_q_1")]

t(apply(d1a[,seq(4,ncol(d1a))], 2, function(x) table(factor(x, levels = c(0,1,NA)), exclude = NULL)))




#                                      0     1   <NA>
#aneurysm_diagnosis_adu_q_1       146769   435   3909
#aneurysm_presence_adu_q_2             0     0 151113
#angioplasty_bypass_adu_q_1       145005  1978   4130
#arrhythmia_presence_adu_q_1      115580 35425    108
#atherosclerosis_presence_adu_q_1      0   701 150412
#carotid_stenosis_adu_q_1         145356   395   5362
#heartattack_presence_adu_q_1     149037  1642    434
#heartfailure_presence_adu_q_1    129253  1110  20750
#heartvalve_presence_adu_q_1           0  1557 149556
#highcholesterol_presence_adu_q_1 110574 19278  21261
#stenosis_presence_adu_q_1             0     0 151113
#stroke_presence_adu_q_1          149214  1175    724



d1a <- d1a[,c("PROJECT_PSEUDO_ID","AGE","GENDER","aneurysm_diagnosis_adu_q_1","angioplasty_bypass_adu_q_1","arrhythmia_presence_adu_q_1", "atherosclerosis_presence_adu_q_1", "carotid_stenosis_adu_q_1","heartattack_presence_adu_q_1","heartfailure_presence_adu_q_1", "heartvalve_presence_adu_q_1", "highcholesterol_presence_adu_q_1","stroke_presence_adu_q_1")]
#write.table(d1a, file = "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v2/1a_q_1_cvd.txt", sep = "\t", quote = F, row.names = F)

d1a <- d1a[,c("PROJECT_PSEUDO_ID","AGE","GENDER","aneurysm_diagnosis_adu_q_1","angioplasty_bypass_adu_q_1","arrhythmia_presence_adu_q_1",  "carotid_stenosis_adu_q_1","heartattack_presence_adu_q_1","heartfailure_presence_adu_q_1",  "highcholesterol_presence_adu_q_1","stroke_presence_adu_q_1")]

healthy = d1a[which(apply(d1a[,seq(4,ncol(d1a))],1,function(x) all(x == 0))),"PROJECT_PSEUDO_ID"]

d <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/1b_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")
d[d == "$6"] <- NA
d[d == "$7"] <- NA
d[d == "3"] <- NA
d[d == ""] <- NA
d[d == "2"] <- 0
d1b <- d[,c("PROJECT_PSEUDO_ID","claudication_followup_adu_q_1","heartattack_followup_adu_q_1", "heartfailure_followup_adu_q_1")]
#write.table(d1b, file = "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v2/1b_q_1_cvd.txt", sep = "\t", quote = F, row.names = F)
t(apply(d1b[,seq(2,ncol(d1b))], 2, function(x) table(factor(x), exclude = NULL)))
new_cvd_1b <- d1b[which(apply(d1b[,seq(2,ncol(d1b))],1,function(x) any(x == 1))),"PROJECT_PSEUDO_ID"]

#                                 0   1   <NA>
#claudication_followup_adu_q_1 6467 548 118170
#heartattack_followup_adu_q_1  6477 382 118326
#heartfailure_followup_adu_q_1 6467 868 117850


d <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/1c_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")
d[d == "$6"] <- NA
d[d == "$7"] <- NA
d[d == "3"] <- NA
d[d == ""] <- NA
d[d == "2"] <- 0
d1c <- d[,c("PROJECT_PSEUDO_ID","claudication_followup_adu_q_1","heartattack_followup_adu_q_1", "heartfailure_followup_adu_q_1","stroke_followup_adu_q_1")]
#write.table(d1c, file = "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v2/1c_q_1_cvd.txt", sep = "\t", quote = F, row.names = F)
t(apply(d1c[,seq(2,ncol(d1c))], 2, function(x) table(factor(x), exclude = NULL)))
new_cvd_1c <- d1c[which(apply(d1c[,seq(2,ncol(d1c))],1,function(x) any(x == 1))),"PROJECT_PSEUDO_ID"]

#                                  0   1  <NA>
#claudication_followup_adu_q_1 63265 617 32216
#heartattack_followup_adu_q_1  63469 341 32288
#heartfailure_followup_adu_q_1 63159 854 32085
#stroke_followup_adu_q_1       63494 279 32325



d <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/2a_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")
d[d == "$6"] <- NA
d[d == "$7"] <- NA
d[d == "3"] <- NA
d[d == ""] <- NA
d[d == "2"] <- 0
d2 <- d[,c("PROJECT_PSEUDO_ID","claudication_followup_adu_q_1","heartattack_followup_adu_q_1", "heartfailure_followup_adu_q_1","stroke_followup_adu_q_1","cvd_followup_adu_q_1", "t2d_followup_adu_q_1")]
write.table(d2, file = "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v2/2a_q_1_cvd.txt", sep = "\t", quote = F, row.names = F)
t(apply(d2[,seq(2,ncol(d2))], 2, function(x) table(factor(x), exclude = NULL)))
new_cvd_2a <- d2[which(apply(d2[,seq(2,ncol(d2))],1,function(x) any(x == 1))),"PROJECT_PSEUDO_ID"]

#                                  0    1  <NA>
#claudication_followup_adu_q_1  9311  498 91576
#heartattack_followup_adu_q_1   9638  448 91299
#heartfailure_followup_adu_q_1  9071  814 91500
#stroke_followup_adu_q_1        9464  291 91630
#cvd_followup_adu_q_1          95404 2300  3681


d <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/3a_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")
d[d == "$6"] <- NA
d[d == "$7"] <- NA
d[d == "3"] <- NA
d[d == ""] <- NA
d[d == "2"] <- 0
d3a <- d[,c("PROJECT_PSEUDO_ID","claudication_followup_adu_q_1","heartattack_followup_adu_q_1", "heartfailure_followup_adu_q_1","stroke_followup_adu_q_1","cvd_followup_adu_q_1", "t2d_followup_adu_q_1")]
#write.table(d3a, file = "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v2/3a_q_1_cvd.txt", sep = "\t", quote = F, row.names = F)
t(apply(d3a[,seq(2,ncol(d3a))], 2, function(x) table(factor(x), exclude = NULL)))
new_cvd_3a <- d3a[which(apply(d3a[,seq(2,ncol(d3a))],1,function(x) any(x == 1))),"PROJECT_PSEUDO_ID"]
 
#                                  0   1  <NA>
#claudication_followup_adu_q_1   627  69 12478
#heartattack_followup_adu_q_1    579 116 12479
#heartfailure_followup_adu_q_1   520 174 12480
#stroke_followup_adu_q_1         631  62 12481
#cvd_followup_adu_q_1          12353 696   125

