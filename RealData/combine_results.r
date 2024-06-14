library(kableExtra)
#combine results
single_cell <- read.csv("single_cell/distance_study.csv")
bbc <- read.csv("bbc/distance_study.csv")
sources <- read.csv("3sources/distance_study.csv")

#table with f score results for different distances used with bisilhouette score
all_dis <- t(rbind(sources, bbc, single_cell)[,2:13])
all_dis <- cbind(rep(c("Euclidean", "Cosine", "Manhattan"), each=4),
            rep(c("E", "C", "M", "F"), 3),
             all_dis)
all_dis <- data.frame(all_dis)
all_dis[,3:8] <- apply(all_dis[,3:8],2, as.numeric)
colnames(all_dis) <- c("","", rep(c("3sources","BBCSport","A549"), each=2))
write.csv(all_dis, "all_results_dis.csv", row.names=FALSE)

res_text <- kbl(all_dis,booktabs=T,"latex",
        escape = FALSE, digits = 4,  
        format.args = list(scientific = FALSE),
        row.names=FALSE)
df <- kable_styling(res_text)
sink("all_results_dis.txt")
print(df)
sink()







# read.csv("ResNMTF/ResNMTF/RealData/bbcsport/distance_study.csv")
# # all <- data.frame("Method" = single_cell$Method, "3Sources" = sources$F.score, "BBCSport" = bbc$F.score, "A549" = single_cell$F.score)

# write.csv(all, "ResNMTF/RealData/all_results_fscores.csv",row.names=FALSE)
# res_text <- kbl(all,booktabs=T,"latex",escape=FALSE,digits=4,row.names=FALSE)
# df <- kable_styling(res_text)
# #save
# sink("ResNMTF/RealData/all_results_fscores.txt")
# print(df)
# sink()
# read.csv("ResNMTF/ResNMTF/RealData/bbcsport/distance_study.csv")
