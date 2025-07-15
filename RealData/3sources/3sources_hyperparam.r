# 3sources analysis
args <- commandArgs(trailingOnly = TRUE)
phi <- as.numeric(args[1])
source("SimStudy/Functions/evaluation_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
library(R.matlab)
library(resnmtf)
library(doParallel)
library(foreach)
library(purrr)
# load in
three_data <- import_matrix("RealData/3sources", "3sources_all_diff")
docs_labs <- read.csv("RealData/3sources/true_labels.csv", row.names = 1)

phi_mat <- matrix(0, 3, 3)
phi_mat[1, c(2, 3)] <- 1
phi_mat[2, 3] <- 1
three_dt <- lapply(three_data, t)

# find optimal values
n_psi <- length(seq(0, 2000, 50))
psi_vec <- rep(seq(0, 2000, 50), 3)
dis_vec <- c(rep("euclidean", n_psi),
             rep("cosine", n_psi),
             rep("manhattan", n_psi))
psi_val <- psi_vec[phi]
dis <- dis_vec[phi]
n_reps <- 5
n_views <- length(three_dt)
n_col <- 3 + 6 * (n_views + 1)

safe_resnmtf <- purrr::possibly(resnmtf::apply_resnmtf, otherwise = NULL)
doParallel::registerDoParallel(min(parallel::detectCores(), 5))
res_list <- foreach::foreach(j = 1:5) %dopar%{
  set.seed(phi + j)
  res_euc <- safe_resnmtf(
    three_data,
    k_min = 4,
    k_max = 8,
    psi = (psi_val) * phi_mat,
    distance = dis,
    stability = FALSE
  )
  if (!is.null(res_euc)) {
    c(
      j, psi_val,
      dis_results(
        three_data, docs_labs, res_euc, psi_val, j,
        paste0("RealData/3sources/", dis)
      )
    )
  } else {
    c(j, psi_val, rep(0, n_col - 2))
  } 
}
results <- rbind(res_list[[1]], res_list[[2]], res_list[[3]], res_list[[4]], res_list[[5]])
write.csv(results, paste0("RealData/3sources/data/three_s_psi_", dis, psi_val, ".csv"))