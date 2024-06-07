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
bbc_rows <- read.csv("bbc/bbc_rows_truth.csv")[,2:6]
bbc_d2 <- import_matrix("bbc/bbc_data_processed.xlsx")
set.seed(10+phi)
phi_bbc <- matrix(0, 2, 2)
phi_bbc[1,2] <- 1
phi_vec <- seq(0, 2000, 50)
phi_val <- phi_vec[phi]
n_views <- length(bbc_d2)
n_col <- 3 + 6 * (n_views + 1)
results_euc <- matrix(0, nrow=5, ncol=n_col)

# res_euc <- restMultiNMTF_run(Xinput = bbc_d2, k_min = 4, 
#                                          k_max = 8, psi = phi_vec[phi]*phi_bbc,
#                                           distance = dis, stability=FALSE)
# results_euc <- c(rep, phi_vec[phi],
#                   dis_results(bbc_d2, bbc_rows, res_euc, phi_vec[phi], rep, paste0("bbc/",dis))) 
# write.csv(results_euc, paste0("bbc/bbc_psi_", dis, phi_vec[phi], "_" ,rep,".csv"))


k <- 1
for(j in 1:5){
     res_euc <- restMultiNMTF_run(Xinput = bbc_d2, k_min = 4, 
                                         k_max = 8, psi = phi_val*phi_bbc,
                                          distance = dis, stability=FALSE)
     results_euc[k,] <- c(j, phi_val,
                                dis_results(bbc_d2, bbc_rows, res_euc, phi_val, j, paste0("bbc/",dis))) 
     write.csv(results_euc, paste0("bbc/data/bbc_psi_",dis,phi_val,".csv"))
     k <- k + 1
} 


# res_bbc <- paste0("bbc/bbc_psi_",dis,".csv")
# psi_vec <- seq(0,2000,250)
# n_views <- length(bbc_d2)
# res_bbc <- res_bbc[,3:12]
# old_names <- c("psi",
#                         paste0("F score (V", 1:n_views, ")"),
#                  paste0("Relevance (V", 1:n_views, ")"),  
#                  paste0("Recovery (V", 1:n_views, ")"),
#                         paste0("BiS - E",c(" (V1)"," (V2)")),
#                         paste0("BiS - C",c(" (V1)"," (V2)")),
#                         paste0("BiS - M",c(" (V1)"," (V2)")),"k")
# colnames(res_bbc) <- old_names

# #calculate means
# measures <- c("F score", "BiS", "Recovery", "Relevance")
# for(measure in measures){
#         m <- paste0(measure,c(" (V1)"," (V2)")) 
#         res_bbc <- cbind(res_bbc, (res_bbc[,m[1]]+res_bbc[,m[2]])/2)
# }

# colnames(res_bbc) <- c(old_names, c("F score", "BiS","Recovery", "Relevance"))

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
# write.csv(final_res, "bbc/final_res_resnmtf.csv")