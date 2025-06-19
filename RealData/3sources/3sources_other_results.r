# local paths
filepath <- "RealData/3sources"
source("SimStudy/OtherMethods/gfa_funcs.r")
source("SimStudy/OtherMethods/biclust_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
source("SimStudy/Functions/evaluation_funcs.r")
library(resnmtf)
three_data <- import_matrix("RealData/3sources", "3sources_all_diff")
docs_labs <- read.csv("RealData/3sources/true_labels.csv", row.names = 1)

# gfa
n_reps <- 5
n_views <- length(three_data)
bbc_res2 <- vector("list", length = n_reps * 2)
n_col <- (n_views + 1) * 6 + 2
results <- matrix(0, nrow = n_reps * 2, ncol = n_col)
col_names <- c(
  "rep", "method",
  paste0("F score (V", 1:n_views, ")"), "F score",
  paste0("Relevance (V", 1:n_views, ")"), "Relevance",
  paste0("Recovery (V", 1:n_views, ")"), "Recovery",
  paste0("BiS-E (V", 1:n_views, ")"), "BiS-E",
  paste0("BiS-C (V", 1:n_views, ")"), "BiS-C",
  paste0("BiS-M (V", 1:n_views, ")"), "BiS-M"
)
colnames(results) <- col_names
# gfa
k <- 1
for (i in 1:n_reps) {
  set.seed(10 + i)
  bbc_res2[[k]] <- gfa_apply(lapply(three_data, t), dim(docs_labs)[2])
  # assess performance
  bisil <- calc_all_sils(three_data, list(
    "row_clusters" = bbc_res2[[k]]$col_clusters,
    "col_clusters" = bbc_res2[[k]]$row_clusters
  ))
  jaccs <- all_jaccs(docs_labs, bbc_res2[[k]]$row_clusters)
  results[k, ] <- c(i, "gfa", jaccs, bisils$euc, bisils$cos, bisils$man)
  k <- k + 1
  write.csv(results, paste0(filepath, "/gfa_method_results.csv"))
}

# biclust methods
set.seed(40)
bcplaid_results <- biclust_results(three_data, docs_labs,
                                   method = BCPlaid, name = "Plaid",
                                   transposed = TRUE)
bcspectral_results <- biclust_results(three_data, docs_labs,
                                      method = BCSpectral, name = "spectral",
                                      transposed = TRUE)
full_res <- rbind(bcplaid_results,
                  bcspectral_results)
colnames(full_res) <- col_names
colnames(bcplaid_results) <- col_names
colnames(bcspectral_results) <- col_names
# save biclust results
write.csv(full_res, "RealData/3sources/biclust_results.csv",
          row.names = FALSE)


results_issvd <- matrix(0, nrow = n_reps, ncol = n_col)
colnames(results_issvd) <- c(
  "rep", "method",
  paste0("F score (V", 1:n_views, ")"), "F score",
  paste0("Relevance (V", 1:n_views, ")"), "Relevance",
  paste0("Recovery (V", 1:n_views, ")"), "Recovery",
  paste0("BiS-E (V", 1:n_views, ")"), "BiS-E",
  paste0("BiS-C (V", 1:n_views, ")"), "BiS-C",
  paste0("BiS-M (V", 1:n_views, ")"), "BiS-M"
)
k <- 1
# save results from python doc
for (t in 1:n_reps) {
    row_clusts <- import_matrix(
    paste0(filepath, "/issvd_res/"),
    paste0(t - 1, "_row_clusts")
  )
  col_clusts <- import_matrix(
    paste0(filepath, "/issvd_res/"),
    paste0(t - 1, "_col_clusts")
  )
  set.seed(10 + t)
  bisils <- calc_all_sils(
    three_data,
    list("row_clusters" = col_clusts, "col_clusters" = row_clusts)
  )
  jaccs <- all_jaccs(docs_labs, row_clusts)
  results_issvd[k, ] <- c(t, "issvd", jaccs, bisils$euc, bisils$cos, bisils$man)
  k <- k + 1
}
write.csv(results_issvd, paste0(filepath, "/python_method_results.csv"))

# mvclustering
results_mvc <- matrix(0, nrow = n_reps, ncol = n_col)
colnames(results_mvc) <- col_names
k <- 1
col_clusts <- list(
  matrix(1, nrow = dim(three_data[[1]])[1], ncol = 6),
  matrix(1, nrow = dim(three_data[[2]])[1], ncol = 6),
  matrix(1, nrow = dim(three_data[[3]])[1], ncol = 6)
)
# save results from python doc
for (t in 1:n_reps) {
  set.seed(10 + t)
  row_clusts <- import_matrix(
    paste0(filepath, "/mvlearn_res/"),
    paste0(t - 1, "_row_clusts")
  )
  bisils <- calc_all_sils(
    three_data,
    list("row_clusters" = col_clusts, "col_clusters" = row_clusts)
  )
  jaccs <- all_jaccs(docs_labs, row_clusts)
  results_mvc[k, ] <- c(t, "issvd", jaccs, bisils$euc, bisils$cos, bisils$man)
  k <- k + 1
}
write.csv(results_mvc, paste0(filepath, "/mvc_method_results.csv"))


# combine best results
choose_best <- function(res, col_name) {
  max_val <- max(res[, col_name])
  res[res[, col_name] == max_val, ]
}
best_gfa <- choose_best(results, "F score")
best_issvd <- choose_best(results_issvd, "F score")
best_mvc <- choose_best(results_mvc, "F score")
best_bcplaid <- choose_best(bcplaid_results, "F score")
best_bcspectral <- choose_best(bcspectral_results, "F score")

# combine results
all_results <- rbind(
  best_gfa,
  best_issvd,
  best_mvc,
  best_bcplaid,
  best_bcspectral
) 
write.csv(all_results, paste0(filepath, "/other_results.csv"), row.names = FALSE)