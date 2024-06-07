args = commandArgs(trailingOnly = TRUE)
dis = as.character(args[1])
n_views = as.numeric(args[2])
phi_vec <- seq(0, 2000, 50)
file_path <- paste0(dis, "/data/", dis, "_psi_")

n_col <- 3 + 6 * (n_views + 1)

results_euc <- matrix(0, nrow=5*length(phi_vec), ncol=n_col)
results_mann <- matrix(0, nrow=5*length(phi_vec), ncol=n_col)
results_cos <- matrix(0, nrow=5*length(phi_vec), ncol=n_col)

k <- 1
for(phi_val in phi_vec){
    file_cosine <- paste0(file_path,"cosine",phi_val,".csv")
    file_euc <- paste0(file_path,"euclidean",phi_val,".csv")
    file_mann <- paste0(file_path,"manhattan",phi_val,".csv")
    results_cos[k:(k+4),] <- as.matrix(read.csv(file_cosine, row.names=1))
    results_euc[k:(k+4),] <- as.matrix(read.csv(file_euc, row.names=1))
    results_mann[k:(k+4),] <- as.matrix(read.csv(file_mann, row.names=1))
    k <- k + 5
}

write.csv(results_euc, paste0(dis, "/data/euc_results.csv"))
write.csv(results_cos, paste0(dis, "/data/cos_results.csv"))
write.csv(results_mann, paste0(dis, "/data/mann_results.csv"))