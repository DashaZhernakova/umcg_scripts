a1 <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/1a_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")
b1 <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/1b_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")
c1 <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/1c_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")
a2 <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/2a_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")
a3 <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/3a_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")
b3 <- read.table("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/3b_q_1_results.csv", header = T, as.is = T, check.names = F, sep = ",")

a1[a1 == "$6"] <- -9
a1[a1 == "3"] <- -9
a1[a1 == ""] <- NA
a1[a1 == "2"] <- 0

#a1$ms_presence_adu_q_1 <- ifelse(is.na(a1$ms_presence_adu_q_1), a1$ams_presence_adu_q_2, a1$ms_presence_adu_q_1)

diseases <- c("ibs", "crohns", "colitisulcerosa", "arthritis", "psoriasis", "parkinsons", "coeliac")
res <- data.frame(matrix(nrow = length(diseases), ncol = 7))
row.names(res) <- diseases
colnames(res) <- c("1a_cases", "1a_missing", "1a_not_asked", "followup_cases", "of_followup_cases_1a_cases", "of_followup_cases_1a_missing", "of_followup_cases_not_asked")
#d = "ibs"
for (d in diseases){

    if(d != "psoriasis"){
        b1_ids <- b1[b1[,paste0(d, "_followup_adu_q_1")] == 1,"project_pseudo_id"]
        c1_ids <- c1[c1[,paste0(d, "_followup_adu_q_1")] == 1,"project_pseudo_id"]
        a2_ids <- a2[a2[,paste0(d, "_followup_adu_q_1")] == 1,"project_pseudo_id"]
        a3_ids <- a3[a3[,paste0(d, "_followup_adu_q_1")] == 1,"project_pseudo_id"]
        b3_ids <- b3[b3[,paste0(d, "_followup_adu_q_1")] == 1,"project_pseudo_id"]

        #paste(length(b1_ids), length(c1_ids), length(a2_ids), length(a3_ids), length(b3_ids))

        followup_positive <- union(union(union(b1_ids, c1_ids), union(a2_ids, a3_ids)), b3_ids)
        #paste(length(followup_positive))
    } else {
        followup_positive <- a3[a3[,"psoriasis_diagnosis_adu_q_1"] == 1,"project_pseudo_id"]
    }

    a1_disease <- a1[,c("project_pseudo_id", paste0(d, "_presence_adu_q_1"))]

    cat("Disease: ", d, "\n")
    cat("Case/controls at baseline:\n")
    print(table(a1_disease[,2], useNA = "always"))
    print(length(followup_positive))
    print(table(a1_disease[a1_disease$project_pseudo_id %in% followup_positive,2], useNA = "always"))
    res[d, 1:3] <- table(a1_disease[,2], useNA = "always")

    res[d,4] <- length(followup_positive)

    res[d,5:7] <- table(a1_disease[a1_disease$project_pseudo_id %in% followup_positive,2], useNA = "always")
    
}