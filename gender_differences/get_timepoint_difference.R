library(RColorBrewer)
library(ggplot2)

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

setwd("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/tmp/")

d <- read.delim("/groups/umcg-lifelines/tmp01/projects/phenotypes/tab_seperated_labels/Laboratory_assessment_Blood.dat", header = T, sep = "\t", as.is = T, check.names = F)
d <- d[d$NUCHTER == "Yes", ]
a1 <- d[d$ENCOUNTERCODE == "Baseline assessment (1A)", ]
a2 <- d[d$ENCOUNTERCODE == "Second assessment (2A)", ]
rownames(a1) <- a1[, 1]
rownames(a2) <- a2[, 1]

a1 <- subset(a1, select = -c(PSEUDOIDEXT, ENCOUNTERCODE, NUCHTER))
a2 <- subset(a2, select = -c(PSEUDOIDEXT, ENCOUNTERCODE, NUCHTER))

a1m <- a1[match(row.names(a2), row.names(a1), nomatch = 0), ]
a2m <- a2[match(row.names(a1m), row.names(a2), nomatch = 0), ]

dif <- a1m - a2m
dif <- dif[,which(sapply(dif, function(x) all(is.na(x))) == F)]

pheno <- read.delim("/groups/umcg-lifelines/tmp01/projects/phenotypes/tab_seperated_labels/Participant (Pa_99_G).dat", header = T, sep = "\t", as.is = T, check.names = F)
rownames(pheno) <- pheno$PSEUDOIDEXT
pheno_m <- pheno[match(rownames(a1m), rownames(pheno), nomatch = 0), ]
all(rownames(pheno_m) == rownames(a1m))
pheno_m <- pheno_m[, c("GESLACHT", "AGE_1A1", "AGE_2A1")]
pheno_m$mean_age <- (pheno_m$AGE_1A1 + pheno_m$AGE_2A1)/2
pheno_m$GESLACHT <- as.factor(pheno_m$GESLACHT)
m <- cbind(pheno_m, dif)

tmp <- m[,c("mean_age", "GESLACHT", "CHO")]
tmp <- na.omit(tmp)
tmp$CHOpos <- ifelse(tmp$CHO > 0, tmp$CHO, NA)
tmp$CHOneg <- ifelse(tmp$CHO < 0, tmp$CHO, NA)
tmp$CHOabs <- abs(tmp$CHO)

male <- aggregate(tmp[tmp$GESLACHT == "Male","CHO"], list(tmp[tmp$GESLACHT == "Male","mean_age"]), mean)
fem <- aggregate(tmp[tmp$GESLACHT == "Female","CHO"], list(tmp[tmp$GESLACHT == "Female","mean_age"]), mean)


pdf("CHO_hex.pdf", width = 5, height = 10)
par(mfrow=c(2,1))
ggplot(tmp[tmp$GESLACHT == "Male",], aes(x = mean_age, y = CHOabs)) + 
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density") +
  geom_point(shape = '.')
ggplot(tmp[tmp$GESLACHT == "Female",], aes(x = mean_age, y = CHOabs)) + 
  stat_density_2d(aes(fill = stat(level)), geom = 'polygon') +
  scale_fill_viridis_c(name = "density") +
  geom_point(shape = '.')
dev.off()

png("CHO_dif_scatter.png",width = 5, height = 5, units = 'in', res = 400)
palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
plot(pheno_m$mean_age, dif$CHO, col = as.factor(pheno_m$GESLACHT))
dev.off()




pheno_w$menopause1 <- NA
pheno_w$menopause2 <- NA
pheno_w[pheno_w$FEM4 == Yes, ] <- 0 # is pregnant
pheno_w[pheno_w$FEM8B == Yes, ] <- NA # hormonal contraceptives last month
pheno_w[pheno_w$FEM7A1 == Yes || pheno_w$FEM7B1 == Yes, pheno_w$FEM7C1 == Yes] <- NA

FEM6C # age of last menstruation


pheno_w$FEM6