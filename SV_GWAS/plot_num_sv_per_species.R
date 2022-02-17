library(ggplot2)
dsv <- read.delim("dSV_per_cohort.strict_filtering.txt", sep = "\t", header = T, as.is = T)
vsv <- read.delim("vSV_per_cohort.txt", sep = "\t", header = T, as.is = T)
dsv$sp <- gsub("\\:[0-9]+", "", dsv[,1])
vsv$sp <- gsub("\\:[0-9]+", "", vsv[,1])
dsv_c <- as.data.frame(table(dsv$sp))
row.names(dsv_c) <- dsv_c[,1]
vsv_c <- as.data.frame(table(vsv$sp))
row.names(vsv_c) <- vsv_c[,1]
counts <- merge(dsv_c, vsv_c, by='row.names', all = T)

counts <- subset(counts, select = c(Row.names,Freq.x, Freq.y))
colnames(counts) <- c("species", "dSV", "vSV")
counts[is.na(counts$dSV), "dSV"] <- 0
counts[is.na(counts$vSV), "vSV"] <- 0
counts$species <- factor(counts$species, levels = counts[order(counts$dSV + counts$vSV, decreasing = T),"species"])
counts2 <- melt(counts)


pdf("num_sv_per_species_barplot.pdf", height = 5, width = 12)
ggplot(counts2, aes(x = species, y = value, fill = variable)) +
    geom_bar(stat="identity", color="grey", position=position_dodge()) +
    scale_fill_brewer(palette="Blues") +
    labs(title="Number of SVs per species after filtering", x = "Species", y = "Number of SVs") +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none")
dev.off()