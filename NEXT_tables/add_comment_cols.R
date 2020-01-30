t <- read.delim("/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/filesreadyforcommentscript/Registration_LL_NEXT_samples_NEXTnumberAdded_PartMerge_v4.0.txt", sep = "\t", quote = "",header = F, as.is = T, check.names = F)
t <- as.data.frame(t)
cs <- which(t[1,] == "Lifelines + NEXT number")
cs <- cs[seq(2,length(cs))] #don't need a comment after row description

empty_col <- rep("", nrow(t))
empty_col[4] = "comment"

res <- t[0:(cs[1] - 1)]
st = cs[1]
for (c in c(cs[2:length(cs)], ncol(t) + 1)){
  res <- cbind(res, empty_col, t[,st:(c - 1)])
  st = c
}
res <- cbind(res, empty_col) # comment for the last LL id
write.table(res, file = "/Users/dashazhernakova/Documents/UMCG/data/NEXT_tables/Registration_LL_NEXT_samples_edSJ_NEXTnumbersAdded_comments.txt", sep = "\t", quote = F, row.names = F, col.names=F)
