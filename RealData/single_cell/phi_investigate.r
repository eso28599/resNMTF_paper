args = commandArgs(trailingOnly = TRUE)
dis = as.character(args[1])
phi = as.numeric(args[2])
#single cell psi selection
source("evaluation_funcs.r")
source("extra_funcs.r")
source("main.r")
source("visualisation.r")
source("gfa_funcs.r")
library(R.matlab)
library(viridis)
library(ggplot2)
library(latex2exp)

labs <- read.csv("single_cell/true_labs.csv")[,2:4]
cell_data <- import_matrix("single_cell/data_processed.xlsx")
phi_mat <- matrix(0,2,2)
phi_mat[1,2] <- 1

n_views <- length(cell_data)
phi_vec <- seq(0, 40000, 500)
phi_val <- phi_vec[phi]
n_reps <- 5
n_col <- (n_views+1)*6+3
results <- matrix(0, nrow=n_reps, ncol=n_col)
colnames(results) <- c("rep", "phi",
                 paste0("F score (V", 1:n_views, ")"), "F score", 
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
                 paste0("BiS-E (V", 1:n_views, ")"), "BiS-E", 
                 paste0("BiS-M (V", 1:n_views, ")"), "BiS-M",
                 paste0("BiS-C (V", 1:n_views, ")"), "BiS-C","k")
k<-1

for(t in 1:n_reps){
    res <- restMultiNMTF_run(Xinput = cell_data, k_min = 3, 
                                            k_max = 6, phi = phi_val*phi_mat, 
                                            distance = dis, stability=FALSE)
    results[k, ] <- c(t, phi_val,
                                dis_results(cell_data, labs, res, phi_val, t, paste0("single_cell/",dis),row_same=TRUE)) 
    k <- k+1
    write.csv(results, paste0("single_cell/data/sc_phi_",dis,phi_val,".csv"))
}


# k <- 1
# for(j in 1:5){
#      rows <- import_matrix(paste0("single_cell/",dis, "/data/row_clusts", phi_val, "_", j, ".xlsx"))
#      cols_og <- import_matrix(paste0("single_cell/",dis, "/data/col_clusts", phi_val, "_", j, ".xlsx"))
#      cols <- list(matrix(1, nrow=ncol(cols_og[[1]]), ncol=ncol(cols_og[[1]])),
#                 matrix(1, nrow=ncol(cols_og[[2]]), ncol=ncol(cols_og[[2]])))
#      res_euc <- list("row_clusters" =rows, "col_clusters" =cols)
#      results[k,] <- c(j, phi_val,
#                                 dis_results(cell_data, labs, res_euc, phi_val, j, "single_cell/sil", row_same=TRUE)) 
#      write.csv(results, paste0("single_cell/data/sil_psi_",dis,phi_val,".csv"))
#      k <- k + 1
# }
