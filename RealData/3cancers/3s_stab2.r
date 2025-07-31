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
n_reps <- 5

# stability
n_views <- length(bbc_d2)
stab <- seq(0.6, 1, 0.05)
stab_vec <- rep(stab, 3)
method_vec <- c(rep("ResNMTF_F", length(stab)),
                rep("NMTF", length(stab)),
                rep("ResNMTF_BiS", length(stab)))
psi_vec <- c(rep(37000, length(stab)),
             rep(0, length(stab)), 
             rep(38000, length(stab)))
method <- method_vec[as.numeric(index)]
omega <- stab_vec[as.numeric(index)]
psi <- psi_vec[as.numeric(index)]


doParallel::registerDoParallel(min(parallel::detectCores(), 5))
res_list <- foreach::foreach(t = 1:n_reps) %dopar%{
  set.seed(20 + psi + t)
  rows <- import_matrix(
    paste0("RealData/3cancers/data/stability_", method, "/"),
    paste0("row_clusts_", t)
  )
  cols <- import_matrix(
    paste0("RealData/3cancers/data/stability_", method, "/"),
    paste0("col_clusts_", t)
  )
  jacc_mat <- read.csv(
    paste0("RealData/3cancers/data/stability_", method,
     "/jacc_mat_", t, ".csv")
  )
  res <- list(
    row_clusters = rows,
    col_clusters = cols
  )
  for (i in 1:n_views) {
      # set clusters not deemed stable to have 0 members
      res$row_clusters[[i]][, jacc_mat[i, ] < omega] <- 0
      res$col_clusters[[i]][, jacc_mat[i, ] < omega] <- 0
  }
  c(t, omega,
      dis_results2(
        bbc_d2, bbc_rows, res, omega, t,
        paste0("RealData/3cancers/", method, "_", t)
      )
  )
}
results <- do.call(rbind, res_list)
print(results)
write.csv(results, paste0("RealData/3cancers/data/stability_", method,
 "/stab_", method, "_", omega, ".csv"))

