library('mgcv')
library('dplyr')
setwd("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/")

pheno_path <- "age_gender_bmi_smk_all_LL.txt"
traits_path <- "pheno_tables/Questionnaire_Health_selected_merged.txt"


col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

merge_match <- function(d1,d2){
        d1_m <- d1[match(row.names(d2), row.names(d1), nomatch = 0 ),]
        d2_m <- d2[match(row.names(d1_m), row.names(d2), nomatch = 0),]
        print(all(row.names(d1_m) == row.names(d2_m)))
        merged <- cbind(d1_m, d2_m)
        colnames(merged) <- c(colnames(d1_m), colnames(d2_m))
        return (merged)
}

# traits of interest
traits0 <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- as.data.frame(sapply(traits0, function(x) as.numeric(as.character(x))))
row.names(traits) <- row.names(traits0)

# Age, gender and other phenotypes
pheno0 <- read.delim(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
pheno <- na.omit(pheno0)

traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))
m <- cbind(pheno_m[,c("age","gender_F1M2")], traits_m)
min_age <- 20
max_age <- 80
m <- m[m$age >= min_age & m$age < max_age,]

prev_w <- aggregate(m[m$gender_F1M2 == 1,], by = list(m[m$gender_F1M2 == 1,"age"]), FUN = mean, na.rm = TRUE)
prev_m <- aggregate(m[m$gender_F1M2 == 2,], by = list(m[m$gender_F1M2 == 2,"age"]), FUN = mean,  na.rm = TRUE)
rownames(prev_w) <- prev_w$age
rownames(prev_m) <- prev_m$age


#prev_w <- aggregate(m[m$gender_F1M2 == 1,], by = list(m[m$gender_F1M2 == 1,"age"]), FUN = sum)
#prev_m <- aggregate(m[m$gender_F1M2 == 2,], by = list(m[m$gender_F1M2 == 2,"age"]), FUN = sum)

# pdf("tmp/diseases_test.pdf", width = 15, height = 15)
# par(mfrow=c(5,4)) 
# for (pheno_idx in seq(4,ncol(traits_m))){
#     tmp_w <- prev_w[, c(1, pheno_idx)]
#     tmp_m <- prev_m[, c(1, pheno_idx)]
#     colnames(tmp_w) <- c("age", "phen")
#     colnames(tmp_m) <- c("age", "phen")
#     plot(x = tmp_w[,1], y = tmp_w[,2], type = 'l', col = "indianred1", main = colnames(traits_m)[pheno_idx])
#     lines(x = tmp_m[,1], y = tmp_m[,2], type = 'l', col = "dodgerblue1")
# }
# dev.off()

#write.table(prev_w, file = "tmp/disease_prevalence_women.txt", sep = "\t", quote = F, col.names = NA)
#write.table(prev_m, file = "tmp/disease_prevalence_men.txt", sep = "\t", quote = F, col.names = NA)


pdf("plots/diseases_prevalence_gam.pdf", width = 15, height = 15)
par(mfrow=c(5,4)) 
for (pheno_idx in seq(4,ncol(prev_m))){
    print(pheno_idx)
    tmp_w <- prev_w[, c(1, pheno_idx)]
    tmp_m <- prev_m[, c(1, pheno_idx)]
    colnames(tmp_w) <- c("age", "phen")
    colnames(tmp_m) <- c("age", "phen")
    gam.w <- gam(phen ~ s(age), data = tmp_w, method = "REML", family = nb)
    gam.m <- gam(phen ~ s(age), data = tmp_m, method = "REML", family = nb)
    pdat.w <- with(tmp_w, expand.grid(age = seq(min_age, max_age, length = 300), 
                                         gender_F1M2 = '1'))
    pdat.w <- transform(pdat.w, pred = predict(gam.w, newdata = pdat.w, type = "response"))
    pdat.m <- with(tmp_m, expand.grid(age = seq(min_age, max_age, length = 300), 
                                         gender_F1M2 = '2'))
    pdat.m <- transform(pdat.m, pred = predict(gam.m, newdata = pdat.m, type = "response"))
    

    merged <- rbind(tmp_w, tmp_m)
    merged$gender_F1M2 <- c(rep(1, nrow(tmp_w)), rep(2, nrow(tmp_m)))
    merged$gender_F1M2 <- as.factor(merged$gender_F1M2)
    merged <- mutate(merged, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
    gam.fit.ordered_sex <- gam(phen ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2), 
                 data = merged, method = 'REML')
    
    gam.p <- summary(gam.fit.ordered_sex)$s.pv[length(summary(gam.fit.ordered_sex)$s.pv)]

    plot(x = tmp_w[,1], y = tmp_w[,2], pch = 21, col = "indianred1", ylim = c(0, max(c(tmp_w$phen, tmp_m$phen), na.rm = T)), xlab = "age", ylab = "disease prevalence", main = paste0(colnames(prev_m)[pheno_idx], "\ninteraction p = ", format(gam.p, digits = 3)))
    lines(pred~age, data = pdat.w, col = col2transparent("indianred1", 80), lwd = 3)

    points(x = tmp_m[,1], y = tmp_m[,2], pch = 21, col = "dodgerblue1")
    lines(pred~age, data = pdat.m, col = col2transparent("dodgerblue1", 80), lwd = 3)
}
dev.off()



pdf("plots/diseases_prevalence_barplot_gam.pdf", width = 15, height = 15)
par(mfrow=c(5,4)) 
for (pheno_idx in seq(4,ncol(prev_m))){
    tmp_w <- prev_w[, c(1, pheno_idx)]
    tmp_m <- prev_m[, c(1, pheno_idx)]
    colnames(tmp_w) <- c("age", "phen")
    colnames(tmp_m) <- c("age", "phen")
    gam.w <- gam(phen ~ s(age), data = tmp_w, method = "REML")
    gam.m <- gam(phen ~ s(age), data = tmp_m, method = "REML")
    pdat.w <- with(tmp_w, expand.grid(age = seq(min_age, max_age - 1), 
                                         gender_F1M2 = '1'))
    pdat.w <- transform(pdat.w, pred = predict(gam.w, newdata = pdat.w, type = "response"))
    pdat.m <- with(tmp_m, expand.grid(age = seq(min_age, max_age - 1), 
                                         gender_F1M2 = '2'))
    pdat.m <- transform(pdat.m, pred = predict(gam.m, newdata = pdat.m, type = "response"))
    

    merged <- rbind(tmp_w, tmp_m)
    merged$gender_F1M2 <- c(rep(1, nrow(tmp_w)), rep(2, nrow(tmp_m)))
    merged$gender_F1M2 <- as.factor(merged$gender_F1M2)
    merged <- mutate(merged, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
    gam.fit.ordered_sex <- gam(phen ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2), 
                 data = merged, method = 'REML')
    
    gam.p <- summary(gam.fit.ordered_sex)$s.pv[length(summary(gam.fit.ordered_sex)$s.pv)]

    rownames(tmp_w) <- tmp_w$age
    rownames(tmp_m) <- tmp_m$age
    x <- merge_match(tmp_w, tmp_m)
    x <- x[c(2,4)]
    colnames(x) <- c("wom", "men")
    x2 <- t(x)
    bp <- barplot(as.matrix(x2), beside = T, border="white", col=c("indianred1","dodgerblue1"),
        xlab = "age", ylab = "disease prevalence", main = paste0(colnames(prev_m)[pheno_idx], "\ninteraction p = ", format(gam.p, digits = 3)))
    mean.bp <- apply(bp,2,mean)
    lines(mean.bp, pdat.w$pred, col = col2transparent("indianred1", 80), lwd = 3)

    lines(mean.bp, pdat.m$pred, col = col2transparent("dodgerblue1", 80), lwd = 3)
}
dev.off()


