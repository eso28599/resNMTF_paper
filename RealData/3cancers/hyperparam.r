# 3sources analysis
args <- commandArgs(trailingOnly = TRUE)
dis <- as.character(args[1])
phi <- as.numeric(args[2])
source("../SimStudy/Functions/evaluation_funcs.r")
source("../SimStudy/Functions/extra_funcs.r")
source("../SimStudy/OtherMethods/gfa_funcs.r")
library(R.matlab)
library(viridis)
library(ggplot2)
library(latex2exp)
library(resnmtf)

# load in
three_data <- import_matrix("3cancers/3sources_all_diff.xlsx")
docs_labs <- read.csv("3cancers/tgca_3cancers_labels.csv", row.names = 1)
set.seed(10 + phi)

phi_mat <- matrix(0, 2, 2)
phi_mat[1, 2] <- 1

# find optimal values
psi_vec <- seq(0, 2000, 50)
psi_vec <- c(0, 200, 500)
psi_val <- psi_vec[phi]
n_reps <- 5
n_views <- length(three_data)
n_col <- 3 + 6 * (n_views + 1)
results <- matrix(0, nrow = n_reps, ncol = n_col)
k <- 1

for (j in 1:1) {
  res_euc <- apply_resnmtf(
    three_data,
    k_min = 3,
    k_max = 6,
    psi = (psi_val) * phi_mat,
    distance = dis,
    stability = FALSE
  )
  results[k, ] <- c(
    j, psi_val,
    dis_results(
      three_data, docs_labs, res_euc, psi_val, j,
      paste0("3sources/", dis)
    )
  )
  write.csv(results, paste0("3cancers/data/three_psi_", dis, psi_val, ".csv"))
  k <- k + 1
}
