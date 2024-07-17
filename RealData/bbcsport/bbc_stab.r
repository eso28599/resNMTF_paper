args = commandArgs(trailingOnly = TRUE)
method = as.character(args[1])
psi = as.numeric(args[2])
#3sources analysis 
source("evaluation_funcs.r")
source("extra_funcs.r")
source("resNMTF_funcs.r")
source("visualisation.r")
source("gfa_funcs.r")
set.seed(10+psi)
bbc_rows <- read.csv("bbc/bbc_rows_truth.csv")[,2:6]
bbc_d2 <- import_matrix("bbc/bbc_data_processed.xlsx")
phi_bbc <- matrix(0, 2, 2)
phi_bbc[1,2] <- 1

#stability
n_views <- length(bbc_d2)
stab_vec <- c("none",seq(0, 1, 0.05))
n_reps <- 5
bbc_res2 <- vector("list", length=n_reps*length(stab_vec))
n_col <- (n_views+1)*6+3
results <- matrix(0, nrow=n_reps*length(stab_vec), ncol=n_col)
colnames(results) <- c("rep", "omega",
                 paste0("F score (V", 1:n_views, ")"), "F score", 
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
                 paste0("BiS-E (V", 1:n_views, ")"), "BiS-E", 
                 paste0("BiS-M (V", 1:n_views, ")"), "BiS-M",
                 paste0("BiS-C (V", 1:n_views, ")"), "BiS-C","k")
k<-1
for(t in 1:n_reps){
        res <- restMultiNMTF_run(Xinput = bbc_d2, k_min = 4, 
                                         k_max = 8, psi = psi*phi_bbc, stab_test=TRUE)
        jacc_mat <- res$jacc
        for(omega in (stab_vec)){
          bbc_res2[[k]] <- res$res
          #perform stability selection
          for (i in 1:n_views){
            #set clusters not deemed stable to have 0 members
            if(omega!="none"){
              bbc_res2[[k]]$row_clusters[[i]][, jacc_mat[i, ] <  omega] <- 0
              bbc_res2[[k]]$col_clusters[[i]][, jacc_mat[i, ] <  omega] <- 0
            }
          }
          results[k,] <-  c(t, omega,
                                dis_results(bbc_d2, bbc_rows, bbc_res2[[k]], omega, t, paste0("bbc/",method))) 
          k <- k+1
        }
        write.csv(results, paste0("bbc/bbc_stab_", method, ".csv"))
}
