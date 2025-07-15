args <- commandArgs(trailingOnly = TRUE)
index <- as.character(args[1])
source("SimStudy/Functions/evaluation_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
library(resnmtf)
library(doParallel)
library(foreach)
library(purrr)

bbc_rows <- read.csv("RealData/single_cell/true_labs.csv", row.names = 1)
bbc_d2 <- import_matrix("RealData/single_cell", "data_processed")
phi_bbc <- matrix(0, 2, 2)
phi_bbc[1,2] <- 1
psi_vec <- c(27000, 0, 32000)
method_vec <- c("ResNMTF_F", "NMTF", "ResNMTF_BiS")
psi <- psi_vec[as.numeric(index)]
method <- method_vec[as.numeric(index)]

# stability
n_views <- length(bbc_d2)
stab_vec <- c("none", seq(0, 1, 0.05))
n_reps <- 5
bbc_res2 <- vector("list", length = n_reps * length(stab_vec))
n_col <- (n_views + 1) * 6 + 3

safe_resnmtf <- purrr::possibly(resnmtf::apply_resnmtf, otherwise = NULL)
doParallel::registerDoParallel(min(parallel::detectCores(), 5))
res_list <- foreach::foreach(t = 1:n_reps) %dopar%{
  set.seed(20 + psi + t)
  res <- safe_resnmtf(
    bbc_d2,
    k_min = 4,
    k_max = 8, phi = psi * phi_bbc, remove_unstable = FALSE, distance = "euclidean", sample_rate = 0.9
  )
  if (is.null(res)) {
    return(c(t, phi_val, rep(0, n_col - 2)))
  } 
  jacc_mat <- res$relevance
  repeat_res <- matrix(0, nrow = length(stab_vec), ncol = n_col)
  k <- 1
  for (omega in (stab_vec)) {
    bbc_res2[[k]] <- res$res
    # perform stability selection
    for (i in 1:n_views) {
      # set clusters not deemed stable to have 0 members
      if (omega != "none") {
        bbc_res2[[k]]$row_clusters[[i]][, jacc_mat[i, ] < omega] <- 0
        bbc_res2[[k]]$col_clusters[[i]][, jacc_mat[i, ] < omega] <- 0
      }
    }
    repeat_res[k, ] <- c(
      t, omega,
      dis_results(
        bbc_d2, bbc_rows, bbc_res2[[k]], omega, t,
        paste0("RealData/3sources/", method), row_same=TRUE
      )
    )
    k <- k + 1
  }
  repeat_res
} 
print(res_list)
results <- do.call(rbind, res_list)
# colnames(results) <- c(
#   "rep", "omega",
#   paste0("F score (V", 1:n_views, ")"), "F score",
#   paste0("Relevance (V", 1:n_views, ")"), "Relevance",
#   paste0("Recovery (V", 1:n_views, ")"), "Recovery",
#   paste0("BiS-E (V", 1:n_views, ")"), "BiS-E",
#   paste0("BiS-M (V", 1:n_views, ")"), "BiS-M",
#   paste0("BiS-C (V", 1:n_views, ")"), "BiS-C", "k"
# )
write.csv(results, paste0("RealData/single_cell/sc_stab_", method, ".csv"))

