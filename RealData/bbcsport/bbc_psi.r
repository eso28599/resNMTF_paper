#3sources analysis 
source("evaluation_funcs.r")
source("extra_funcs.r")
source("resNMTF_funcs.r")
source("visualisation.r")
source("gfa_funcs.r")
library(R.matlab)
bbc_rows <- read.csv("bbc/bbc_rows_truth.csv")[,2:6]
bbc_d2 <- import_matrix("bbc/bbc_data_processed.xlsx")
set.seed(40)
phi_bbc <- matrix(0, 2, 2)
phi_bbc[1,2] <- 1
phi_vec <- c(seq(0,1000,250), 1250, 1500, 1750, 2000)
bbc_res <- vector("list", length=5*length(phi_vec))
results <- matrix(0, nrow=5*length(phi_vec), ncol=11)
k <- 1
for(i in 1:length(phi_vec)){
        for(j in 1:5){
                res <- restMultiNMTF_run(Xinput = bbc_d2, k_min = 4, 
                                         k_max = 8, psi = (phi_vec[i])*phi_bbc, stability=FALSE)
                bbc_res[[k]] <- res
                bisils <- c(sil_score(as.matrix(bbc_d2[[1]]), res$row_clusters[[1]], res$col_clusters[[1]])$sil, 
                        sil_score(as.matrix(bbc_d2[[2]]), res$row_clusters[[2]], res$col_clusters[[2]])$sil)
                jaccs <-  list(jaccard_row(res$col_clusters[[1]],bbc_rows),
                                jaccard_row(res$col_clusters[[2]], bbc_rows))
                results[k,] <- c(i, phi_vec[i], jaccs[[1]]$f_score[1], jaccs[[2]]$f_score[1], 
                                jaccs[[1]]$rel[1], jaccs[[2]]$rel[1], 
                                jaccs[[1]]$rec[1], jaccs[[2]]$rec[1], 
                                bisils, dim(res$row_clusters[[1]])[2])
                write.csv(results, "bbc/bbc_psi.csv")
                k <- k + 1
        }  
}

res_bbc <- read.csv("bbc/bbc_psi.csv")
psi_vec <- seq(0,2000,250)
n_views <- length(bbc_d2)
res_bbc <- res_bbc[,3:12]
old_names <- c("psi",
                        paste0("F score (V", 1:n_views, ")"),
                 paste0("Relevance (V", 1:n_views, ")"),  
                 paste0("Recovery (V", 1:n_views, ")"),
                        paste0("BiS",c(" (V1)"," (V2)")),"k")
colnames(res_bbc) <- old_names

#calculate means
measures <- c("F score", "BiS", "Recovery", "Relevance")
for(measure in measures){
        m <- paste0(measure,c(" (V1)"," (V2)")) 
        res_bbc <- cbind(res_bbc, (res_bbc[,m[1]]+res_bbc[,m[2]])/2)
}

colnames(res_bbc) <- c(old_names, c("F score", "BiS","Recovery", "Relevance"))

#save results for max bis and max f score
max_bis <- which.max(res_bbc[,"BiS"])
max_f <- which.max(res_bbc[,"F score"])
max_bis_res <- res_bbc[max_bis,]
max_f_res <- res_bbc[max_f, "F score"]
nmtf <- which.max(res_bbc[res_bbc[,"psi"]==0,"F score"])

final_res <- res_bbc[c(max_bis,max_f,nmtf),
                c(paste0("F score (V", 1:n_views, ")"), "F score", 
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
                 paste0("BiS (V", 1:n_views, ")"), "BiS")]

final_res <- cbind(c("ResNMTF - BiS", "ResNMTF - F score", "NMTF"), final_res)
colnames(final_res) <- c("Method", paste0("F score (V", 1:n_views, ")"), "F score", 
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
                 paste0("BiS-E (V", 1:n_views, ")"), "BiS-E")
write.csv(final_res, "bbc/final_res_resnmtf.csv")