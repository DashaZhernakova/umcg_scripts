library(RColorBrewer)
library(ggplot2)

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

setwd("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/tmp/")

d <- read.delim("/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v1/tab_separated_labels/Laboratory_assessment_Blood.dat", header = T, sep = "\t", as.is = T, check.names = F)
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

pheno <- read.delim("../age_gender_medications_newCVD_all_LL.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
pheno <- pheno[pheno$statins == 0 && pheno$estrogens ==0 && pheno$antihypertensives == 0,]

ind <- intersect(row.names(dif), row.names(pheno))
pheno_m <- pheno[ind,]
dif_m <- dif[ind,]
all(rownames(pheno_m) == rownames(dif_m))

pheno_m <- subset(pheno_m, select = -c(statins, estrogens, antihypertensives))

m <- cbind(pheno_m, dif_m[,"CHO"])
colnames(m) <- c(colnames(pheno_m), "CHO_dif")



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