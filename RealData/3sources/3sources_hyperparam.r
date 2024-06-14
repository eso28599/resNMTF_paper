#3sources analysis 
args = commandArgs(trailingOnly = TRUE)
dis = as.character(args[1])
phi = as.numeric(args[2])
source("evaluation_funcs.r")
source("extra_funcs.r")
source("main.r")
source("visualisation.r")
source("gfa_funcs.r")
library(R.matlab)
library(viridis)
library(ggplot2)
library(latex2exp)

#load in 
three_data <- import_matrix("3sources/3sources_all_diff.xlsx")
docs_labs <- read.csv("3sources/true_labels.csv", row.names=1)
set.seed(10+phi)

phi_mat <- matrix(0, 3, 3)
phi_mat[1, c(2,3)] <- 1
phi_mat[2, 3] <- 1
three_dt <- lapply(three_data, t)

#find optimal values
psi_vec <- seq(0, 2000, 50)
psi_val <- psi_vec[phi]
n_reps <- 5
n_views <- length(three_dt)
n_col <- 3 + 6 * (n_views + 1)
results <- matrix(0, nrow=n_reps, ncol=n_col)
k <- 1

for(j in 1:5){
     res_euc <- restMultiNMTF_run(Xinput = three_data, k_min = 4, 
                                         k_max = 8,
                                psi = (psi_val)*phi_mat, 
                                distance = dis, 
                                stability = FALSE)
     results[k,] <- c(j, psi_val,
                                dis_results(three_data, docs_labs, res_euc, psi_val, j, paste0("3sources/",dis))) 
     write.csv(results, paste0("3sources/data/three_s_psi_",dis,psi_val,".csv"))
     k <- k + 1
}
