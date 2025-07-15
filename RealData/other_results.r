# local paths
args <- commandArgs(trailingOnly = TRUE)
dataset <- as.character(args[1])
source("SimStudy/OtherMethods/gfa_funcs.r")
source("SimStudy/OtherMethods/biclust_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
source("SimStudy/Functions/evaluation_funcs.r")
library(resnmtf)
library(biclust)
library(foreach)
library(doParallel)
if (dataset == "3sources") {
  filepath <- "RealData/3sources"
  input_data <- import_matrix("RealData/3sources", "3sources_all_diff")
  labels <- read.csv("RealData/3sources/true_labels.csv", row.names = 1)
  transposed_data <- TRUE
} else if (dataset == "3cancers") {
  filepath <- "RealData/3cancers"
  input_data <- import_matrix("RealData/3cancers", "tgca_3cancers")
  input_data <- lapply(input_data, function(x) x[, 2:1242])
  input_data <- lapply(input_data, function(x) apply(x, 2, as.numeric))
  labels <- read.csv("RealData/3cancers/tgca_labels.csv", row.names = NULL)
  transposed_data <- TRUE
} else if (dataset == "bbcsport") {
  filepath <- "RealData/bbcsport"
  labels <- read.csv(paste0(filepath, "/bbc_rows_truth.csv"))[, 2:6]
  input_data <- import_matrix(filepath, "bbc_data_processed")
  transposed_data <- TRUE
} else if (dataset == "single_cell") {
  filepath <- "RealData/single_cell"
  labels <- read.csv(paste0(filepath, "/true_labs.csv"))[, 2:4]
  input_data <- import_matrix("RealData/single_cell", "data_processed")
  transposed_data <- FALSE
} else {
  stop("Invalid dataset specified.")
}

# gfa
n_reps <- 5
n_views <- length(input_data)
bbc_res2 <- vector("list", length = n_reps * 2)
n_col <- (n_views + 1) * 6 + 2
results <- matrix(0, nrow = n_reps * 2, ncol = n_col)
col_names_list <- c(
  "rep", "method",
  paste0("F score (V", 1:n_views, ")"), "F score",
  paste0("Relevance (V", 1:n_views, ")"), "Relevance",
  paste0("Recovery (V", 1:n_views, ")"), "Recovery",
  paste0("BiS-E (V", 1:n_views, ")"), "BiS-E",
  paste0("BiS-C (V", 1:n_views, ")"), "BiS-C",
  paste0("BiS-M (V", 1:n_views, ")"), "BiS-M"
)
colnames(results) <- col_names_list
# # gfa
# doParallel::registerDoParallel(min(parallel::detectCores(), 5))
# res_list <- foreach::foreach(i = 1:n_reps) %dopar% {
res_list <- vector("list", length = n_reps)
for (i in 1:n_reps) {
  set.seed(10 + i)
  if (transposed_data) {
    # if data is transposed, we need to transpose it back for gfa
    gfa_input_data <- lapply(input_data, t)
  } else {
    gfa_input_data <- input_data
  }
  gfa_res <- gfa_apply(gfa_input_data, dim(labels)[2])
  # assess performance
  if (transposed_data) {
    clust_res <- list(
      "row_clusters" = gfa_res$col_clusters,
      "col_clusters" = gfa_res$row_clusters
    )
  } else {
    clust_res <- list(
      "row_clusters" = gfa_res$row_clusters,
      "col_clusters" = gfa_res$col_clusters
    )
  }
  bisils <- calc_all_sils(input_data, clust_res)
  jaccs <- all_jaccs(labels, gfa_res$row_clusters)
  res_list[[i]] <- c(i, "gfa", jaccs, bisils$euc, bisils$cos, bisils$man)
  print(paste0("GFA iteration ", i, " completed."))
}
results <- rbind(
  res_list[[1]], res_list[[2]], res_list[[3]],
  res_list[[4]], res_list[[5]]
)
write.csv(results, paste0(filepath, "/gfa_method_results.csv"))

# biclust methods
set.seed(40)
bcplaid_results <- biclust_results(input_data, labels,
  BCPlaid(),
  name = "Plaid",
  transposed = transposed_data
)
results <- apply_biclust_method(
  input_data, BCPlaid, "Plaid",
  TRUE
)
colnames(bcplaid_results) <- col_names_list
write.csv(bcplaid_results, paste0(filepath, "/bcplaid_results.csv"),
  row.names = FALSE
)
bcspectral_results <- biclust_results(input_data, labels,
  method = BCSpectral, name = "spectral",
  transposed = transposed_data
)
colnames(bcspectral_results) <- col_names_list
# save biclust results
write.csv(bcspectral_results, paste0(filepath, "/bcspectral_results.csv"),
  row.names = FALSE
)

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
  if (transposed_data) {
    clust_res <- list(
      "row_clusters" = col_clusts,
      "col_clusters" = row_clusts
    )
  } else {
    clust_res <- list(
      "row_clusters" = row_clusts,
      "col_clusters" = col_clusts
    )
  }
  bisils <- calc_all_sils(
    input_data,
    clust_res
  )
  jaccs <- all_jaccs(labels, clust_res$row_clusters)
  results_issvd[k, ] <- c(t, "issvd", jaccs, bisils$euc, bisils$cos, bisils$man)
  k <- k + 1
}
write.csv(results_issvd, paste0(filepath, "/python_method_results.csv"))

# mvclustering
results_mvc <- matrix(0, nrow = n_reps, ncol = n_col)
colnames(results_mvc) <- col_names_list
k <- 1
if (length(input_data) != 3) {
  col_clusts <- list(
    matrix(1, nrow = dim(input_data[[1]])[1], ncol = ncol(labels)),
    matrix(1, nrow = dim(input_data[[2]])[1], ncol = ncol(labels))
  )
} else {
  col_clusts <- list(
    matrix(1, nrow = dim(input_data[[1]])[1], ncol = ncol(labels)),
    matrix(1, nrow = dim(input_data[[2]])[1], ncol = ncol(labels)),
    matrix(1, nrow = dim(input_data[[3]])[1], ncol = ncol(labels))
  )
}

# save results from python doc
for (t in 1:n_reps) {
  set.seed(10 + t)
  row_clusts <- import_matrix(
    paste0(filepath, "/mvlearn_res/"),
    paste0(t - 1, "row_mvclusts")
  )
  if (transposed_data) {
    clust_res <- list(
      "row_clusters" = col_clusts,
      "col_clusters" = row_clusts
    )
  } else {
    clust_res <- list(
      "row_clusters" = row_clusts,
      "col_clusters" = col_clusts
    )
  }
  # bisils <- calc_all_sils(
  #   input_data,
  #   clust_res
  # )
  jaccs <- all_jaccs(labels, clust_res$row_clusters)
  results_mvc[k, ] <- c(t, "mvc", jaccs, rep(0, (n_views + 1) * 3))
  k <- k + 1
}
write.csv(results_mvc, paste0(filepath, "/mvc_method_results.csv"))

# combine best results
choose_best <- function(res, col_name, col_names = col_names_list) {
  colnames(res) <- col_names
  max_val <- max(res[, col_name])
  as.matrix(res[res[, col_name] == max_val, ])
}
results <- read.csv(paste0(filepath, "/gfa_method_results.csv"), row.names = 1)
best_gfa <- choose_best(results, "F score")
results_issvd <- read.csv(paste0(filepath, "/python_method_results.csv"),
  row.names = 1
)
best_issvd <- choose_best(results_issvd, "F score")
results_mvc <- read.csv(paste0(filepath, "/mvc_method_results.csv"),
  row.names = 1
)
best_mvc <- choose_best(results_mvc, "F score")
bcplaid_results <- read.csv(paste0(filepath, "/bcplaid_results.csv"))
best_bcplaid <- choose_best(bcplaid_results, "F score")
bcspectral_results <- read.csv(paste0(filepath, "/bcspectral_results.csv"))
best_bcspectral <- choose_best(bcspectral_results, "F score")

# combine results
all_results <- rbind(
  best_gfa,
  best_issvd,
  best_mvc,
  best_bcplaid,
  best_bcspectral
)
write.csv(all_results, paste0(filepath, "/other_results.csv"),
  row.names = FALSE
)
