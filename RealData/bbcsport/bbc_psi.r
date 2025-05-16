# 3sources analysis
args <- commandArgs(trailingOnly = TRUE)
dis <- as.character(args[1])
phi <- as.numeric(args[2])
source("SimStudy/Functions/evaluation_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
source("SimStudy/OtherMethods/gfa_funcs.r")
library(resnmtf)
library(R.matlab)
bbc_rows <- read.csv("RealData/bbcsport/bbc_rows_truth.csv")[, 2:6]
bbc_d2 <- import_matrix("RealData/bbcsport/bbc_data_processed.xlsx")
set.seed(10 + phi)
phi_bbc <- matrix(0, 2, 2)
phi_bbc[1, 2] <- 1
phi_vec <- seq(0, 2000, 50)
phi_val <- phi_vec[phi]
n_views <- length(bbc_d2)
n_col <- 3 + 6 * (n_views + 1)
results_euc <- matrix(0, nrow = 5, ncol = n_col)

k <- 1
for (j in 1:5) {
  res_euc <- restMultiNMTF_run(
    Xinput = bbc_d2, k_min = 4,
    k_max = 8, psi = phi_val * phi_bbc,
    distance = dis, stability = FALSE
  )
  results_euc[k, ] <- c(
    j, phi_val,
    dis_results(bbc_d2, bbc_rows, res_euc, phi_val, j, paste0("bbc/", dis))
  )
  write.csv(results_euc, paste0(
    "RealData/bbcsport/data/bbc_psi_",
    dis, phi_val, ".csv"
  ))
  k <- k + 1
}
