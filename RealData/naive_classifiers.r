source("SimStudy/Functions/evaluation_funcs.r")

# 3sources
n_clusts <- 6
rows <- read.csv("RealData/3sources/true_labels.csv", header = T, row.names = 1)
results <- c()
for (i in 1:100) {
  results <- c(results, jaccard_row(rows, rows[
    sample(169),
    sample(n_clusts)
  ])$f_score[1])
}
three_s <- mean(results)

# A549
n_clusts <- 3
rows <- read.csv("RealData/single_cell/true_labs.csv",
  header = T, row.names = 1
)
results <- c()
for (i in 1:100) {
  results <- c(
    results,
    jaccard_row(rows, rows[sample(2641), sample(n_clusts)])$f_score[1]
  )
}
sc_res <- mean(results)

# BBCSport
n_clusts <- 5
rows <- read.csv("RealData/bbcsport/bbc_rows_truth.csv",
  header = T, row.names = 1
)
results <- c()
for (i in 1:100) {
  results <- c(
    results,
    jaccard_row(rows, rows[sample(544), sample(n_clusts)])$f_score[1]
  )
}
bbcsport <- mean(results)
write.csv(data.frame(
  "3sources" = three_s,
  "A549" = sc_res,
  "BBCSport" = bbcsport
), "RealData/naive_classifiers.csv", row.names = FALSE)
