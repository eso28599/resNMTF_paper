args <- commandArgs(trailingOnly = TRUE)
index <- as.character(args[1])
source("SimStudy/Functions/evaluation_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
library(resnmtf)
library(doParallel)
library(foreach)
library(purrr)

three_data <- import_matrix("RealData/3cancers", "tgca_3cancers")
three_data <- lapply(three_data, function(x) x[, 2:1242])
bbc_d2 <- lapply(three_data, function(x) apply(x, 2, as.numeric))
bbc_rows <- read.csv("RealData/3cancers/tgca_labels.csv", row.names = NULL)
phi_bbc <- matrix(0, 2, 2)
phi_bbc[1,2] <- 1
psi_vec <- c(37000, 0, 38000)
method_vec <- c("ResNMTF_F", "NMTF", "ResNMTF_BiS")
psi <- psi_vec[as.numeric(index)]
method <- method_vec[as.numeric(index)]
stability_path <- paste0("RealData/3cancers/data/stability_", method)

# stability
n_views <- length(bbc_d2)
n_reps <- 5
n_col <- (n_views + 1) * 6 + 3

safe_resnmtf <- purrr::possibly(resnmtf::apply_resnmtf, otherwise = NULL)
doParallel::registerDoParallel(min(parallel::detectCores(), 5))
res_list <- foreach::foreach(t = 1:n_reps) %dopar%{
  set.seed(20 + psi + t)
  res <- safe_resnmtf(
    bbc_d2,
    k_min = 4,
    k_max = 8, psi = psi * phi_bbc, remove_unstable = FALSE, distance = "cosine", sample_rate = 0.9
  )
  jacc_mat <- res$relevance
  results <- res$res
  # save row and column clusters
  export_matrix_list(
    results$row_clusters, paste0("row_clusts_", t), stability_path
  )
  export_matrix_list(
    results$col_clusters, paste0("col_clusts_", t), stability_path
  )
  write.csv(
    jacc_mat, paste0(stability_path, "/jacc_mat_", t, ".csv"), row.names = FALSE
  )
  write.csv(
    results$All_Error, paste0(stability_path, "/errors_", t, ".csv"), row.names = FALSE
  )
} 
