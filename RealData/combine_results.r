library(kableExtra)
# combine results
single_cell <- read.csv("RealData/single_cell/distance_study.csv")
bbc <- read.csv("RealData/bbcsport/distance_study.csv")
sources <- read.csv("RealData/3sources/distance_study.csv")
cancers <- read.csv("RealData/3cancers/distance_study.csv")

# table with f score results for different distances
# used with bisilhouette score
all_dis <- t(rbind(sources, bbc, single_cell, cancers)[, 2:13])
all_dis <- cbind(
  rep(c("Euclidean", "Cosine", "Manhattan"), each = 4),
  rep(c("E", "C", "M", "F"), 3),
  all_dis
)
all_dis <- data.frame(all_dis)
all_dis[, 3:11] <- apply(all_dis[, 3:11], 2, as.numeric)
colnames(all_dis) <- c(
  "", "",
  rep(c("3sources", "BBCSport", "A549", "3cancers"), each = 3)
)
write.csv(all_dis, "RealData/all_results_dis.csv", row.names = FALSE)

res_text <- kbl(all_dis,
  booktabs = T, "latex",
  escape = FALSE, digits = 4,
  format.args = list(scientific = FALSE),
  row.names = FALSE
)
df <- kable_styling(res_text)
sink("RealData/all_results_dis.txt")
print(df)
sink()
