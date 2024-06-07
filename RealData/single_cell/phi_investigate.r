args = commandArgs(trailingOnly = TRUE)
dis = as.character(args[1])
phi = as.numeric(args[2])
#single cell psi selection
source("evaluation_funcs.r")
source("extra_funcs.r")
# source("resNMTF_funcs.r")
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
# res_list <- vector("list", length=n_reps*length(phi_vec))
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
    # res_list[[k]] <- res
            #assess performance 
            # bisils <- calc_all_sils(cell_data, res)
            # jaccs <-  all_jaccs(labs, res$row_clusters)
    # results[k, ] <- c(t, val, jaccs, bisils$euc, bisils$man, bisils$cosine)
    results[k, ] <- c(t, phi_val,
                                dis_results(cell_data, labs, res, phi_val, t, paste0("single_cell/",dis),row_same=TRUE)) 
    k <- k+1
    write.csv(results, paste0("single_cell/data/sc_phi_",dis,phi_val,".csv"))
    # write.csv(results, "single_cell/pi.csv")
}
        


# res_bbc <- read.csv("investigate_applications/single_cell/phi.csv")
# #make dataframe
# psi_vec <- seq(0,2000,250)
# res_bbc <- res_bbc[,3:15]
# n_views <- 2
# old_names <- c("psi",
#                      paste0("F score (V", 1:n_views, ")"), "F score", 
#                  paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
#                  paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
#                  paste0("BiS (V", 1:n_views, ")"), "BiS")
# colnames(res_bbc) <- old_names

# #save results for max bis and max f score
# max_bis <- which.max(res_bbc[,"BiS"])
# max_f <- which.max(res_bbc[,"F score"])
# max_bis_res <- res_bbc[max_bis,]
# max_f_res <- res_bbc[max_f, "F score"]
# nmtf <- which.max(res_bbc[res_bbc[,"psi"]==0,"F score"])

# final_res <- res_bbc[c(max_bis,max_f,nmtf),
#                 c(paste0("F score (V", 1:n_views, ")"), "F score", 
#                  paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
#                  paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
#                  paste0("BiS (V", 1:n_views, ")"), "BiS")]


# final_res <- cbind(c("ResNMTF - BiS", "ResNMTF - F score", "NMTF"), final_res)
# colnames(final_res) <- c("Method", paste0("F score (V", 1:n_views, ")"), "F score", 
#                  paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
#                  paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
#                  paste0("BiS-E (V", 1:n_views, ")"), "BiS-E")
# write.csv(final_res, "investigate_applications/single_cell/final_res_resnmtf.csv")
