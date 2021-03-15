args <- commandArgs(trailingOnly = TRUE)


d1 <- read.delim(args[1], header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
d2 <- read.delim(args[2], header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)


merge_match <- function(d1,d2){
        d1_m <- d1[match(row.names(d2), row.names(d1), nomatch = 0 ),]
        d2_m <- d2[match(row.names(d1_m), row.names(d2), nomatch = 0),]
        print(all(row.names(d1_m) == row.names(d2_m)))
        merged <- cbind(d1_m, d2_m)
        colnames(merged) <- c(colnames(d1_m), colnames(d2_m))
        return (merged)
}

m <- merge_match(d1,d2)
write.table(m, file = args[3], sep = "\t", quote = F, col.names = NA)

