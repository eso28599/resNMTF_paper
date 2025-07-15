args <- commandArgs(trailingOnly = TRUE)
phi <- as.numeric(args[1])
# single cell psi selection
source("SimStudy/Functions/evaluation_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
library(R.matlab)
library(resnmtf)
library(doParallel)
library(foreach)
library(purrr)

labs <- read.csv("RealData/single_cell/true_labs.csv")[, 2:4]
cell_data <- import_matrix("RealData/single_cell/", "data_processed")
phi_mat <- matrix(0, 2, 2)
phi_mat[1, 2] <- 1

n_views <- length(cell_data)
n_psi <- length(seq(0, 40000, 1000))
phi_vec <- rep(seq(0, 40000, 1000), 3)
dis_vec <- c(rep("euclidean", n_psi),
             rep("cosine", n_psi),
             rep("manhattan", n_psi))
phi_val <- phi_vec[phi]
dis <- dis_vec[phi]
n_reps <- 5
n_col <- (n_views + 1) * 6 + 3

safe_resnmtf <- purrr::possibly(resnmtf::apply_resnmtf, otherwise = NULL)
doParallel::registerDoParallel(min(parallel::detectCores(), 5))
res_list <- foreach::foreach(j = 1:5) %dopar%{
    set.seed(10 + phi + j)
    res <- safe_resnmtf(
      cell_data,
      k_min = 3,
      k_max = 6, phi = phi_val * phi_mat,
      distance = dis, stability = FALSE
    )
    if (is.null(res)) {
      return(c(j, phi_val, rep(0, n_col - 2)))
    } else {
      c(
        j, phi_val,
        dis_results(cell_data, labs, res, phi_val, j,
          paste0("single_cell/", dis),
          row_same = TRUE
        )
      )
    }
}

results <- rbind(res_list[[1]], res_list[[2]], res_list[[3]], res_list[[4]], res_list[[5]])
write.csv(results, paste0("RealData/single_cell/data/sc_psi_", dis, phi_val, ".csv"))

