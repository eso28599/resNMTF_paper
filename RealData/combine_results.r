#combine results
single_cell <- read.csv("ResNMTF/RealData/single_cell/all_results.csv")
bbc <- read.csv("ResNMTF/RealData/bbcsport/all_results.csv")
sources <- read.csv("ResNMTF/RealData/3sources/all_results.csv")
all <- data.frame("Method" = single_cell$Method, "3Sources" = sources$F.score, "BBCSport" = bbc$F.score, "A549" = single_cell$F.score)

write.csv(all, "ResNMTF/RealData/all_results_fscores.csv",row.names=FALSE)
res_text <- kbl(all,booktabs=T,"latex",escape=FALSE,digits=4,row.names=FALSE)
df <- kable_styling(res_text)
#save
sink("ResNMTF/RealData/all_results_fscores.txt")
print(df)
sink()
