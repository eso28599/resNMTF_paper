args <- commandArgs(trailingOnly = TRUE)
phi <- as.numeric(args[1])
source("SimStudy/Functions/evaluation_funcs.r")
source("SimStudy/OtherMethods/gfa_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
library(R.matlab)
library(resnmtf)
library(doParallel)
library(foreach)
library(purrr)
# load in
# three_data <- import_matrix("RealData/3cancers/tgca_3cancers.xlsx")
three_data <- import_matrix("RealData/3cancers", "tgca_3cancers")
three_data <- lapply(three_data, function(x) x[, 2:1242])

three_data <- lapply(three_data, function(x) apply(x, 2, as.numeric))
docs_labs <- read.csv("RealData/3cancers/tgca_labels.csv", row.names = NULL)
set.seed(10 + phi)

phi_mat <- matrix(0, 2, 2)
phi_mat[1, 2] <- 1

# find optimal values
n_psi <- length(seq(0, 40000, 1000))
psi_vec <- rep(seq(0, 40000, 1000), 3)
dis_vec <- c(rep("euclidean", n_psi),
             rep("cosine", n_psi),
             rep("manhattan", n_psi))
psi_val <- psi_vec[phi]
dis <- dis_vec[phi]
n_reps <- 5
n_views <- length(three_data)
n_col <- 3 + 6 * (n_views + 1)
results <- matrix(0, nrow = n_reps, ncol = n_col)
k <- 1

safe_resnmtf <- purrr::possibly(resnmtf::apply_resnmtf, otherwise = NULL)
doParallel::registerDoParallel(min(parallel::detectCores(), 5))
res_list <- foreach::foreach(j = 1:5) %dopar%{
  set.seed(10 + phi + j)
  res_euc <- safe_resnmtf(
    three_data,
    k_min = 3,
    k_max = 6,
    psi = (psi_val) * phi_mat,
    distance = dis,
    stability = FALSE
  )
  if (!is.null(res_euc)) {
    c(
    j, psi_val,
    dis_results(
      three_data, docs_labs, res_euc, psi_val, j,
      "RealData/3sources"
    )
  )
  } else {
    c(j, psi_val, rep(0, n_col - 2))
  }
}
results <- rbind(res_list[[1]], res_list[[2]], res_list[[3]], res_list[[4]], res_list[[5]])
write.csv(results, paste0("RealData/3cancers/data/three_psi_", dis, psi_val, ".csv"))