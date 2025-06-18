args <- commandArgs(trailingOnly = TRUE)
phi <- as.numeric(args[1])
source("SimStudy/Functions/evaluation_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
source("SimStudy/OtherMethods/gfa_funcs.r")
library(resnmtf)
library(R.matlab)
library(doParallel)
library(foreach)
bbc_rows <- read.csv("RealData/bbcsport/bbc_rows_truth.csv")[, 2:6]
bbc_d2 <- import_matrix("RealData/bbcsport", "bbc_data_processed")
set.seed(10 + phi)
phi_bbc <- matrix(0, 2, 2)
phi_bbc[1, 2] <- 1
n_psi <- length(seq(0, 2000, 50))
psi_vec <- rep(seq(0, 2000, 50), 3)
dis_vec <- c(rep("euclidean", n_psi),
             rep("cosine", n_psi),
             rep("manhattan", n_psi))
phi_val <- psi_vec[phi]
dis <- dis_vec[phi]
n_views <- length(bbc_d2)
n_col <- 3 + 6 * (n_views + 1)
results_euc <- matrix(0, nrow = 5, ncol = n_col)

doParallel::registerDoParallel(min(parallel::detectCores(), 5))
res_list <- foreach::foreach(j = 1:5) %dopar%{
  res_euc <- apply_resnmtf(
    bbc_d2, k_min = 4,
    k_max = 8, psi = phi_val * phi_bbc,
    distance = dis, stability = FALSE
  )
  c(
    j, phi_val,
    dis_results(bbc_d2, bbc_rows, res_euc, phi_val, j, paste0("RealData/bbcsport/", dis))
  )
}
print(res_list)
results <- rbind(res_list[[1]], res_list[[2]], res_list[[3]],
                 res_list[[4]], res_list[[5]])
write.csv(results_euc, paste0(
    "RealData/bbcsport/data/bbc_psi_",
    dis, phi_val, ".csv"
))
