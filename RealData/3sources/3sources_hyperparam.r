#3sources analysis 
#sim_common_files/resNMTF_files/
source("evaluation_funcs.r")
source("extra_funcs.r")
source("resNMTF_funcs.r")
source("visualisation.r")
source("gfa_funcs.r")
library(R.matlab)
library(viridis)
library(ggplot2)
library(latex2exp)

#load in 
three_data <- import_matrix("3sources/3sources_all_diff.xlsx")
docs_labs <- read.csv("3sources/true_labels.csv", row.names=1)
set.seed(10)

phi_mat <- matrix(0, 3, 3)
phi_mat[1, c(2,3)] <- 1
phi_mat[2, 3] <- 1
three_dt <- lapply(three_data, t)

#find optimal values
psi_vec <- c(seq(0,1000,250), 1250, 1500, 1750, 2000)
n_reps <- 5
three_run <- vector("list", length=length(psi_vec)*n_reps)
results <- matrix(0, nrow=n_reps*length(psi_vec), ncol=17)
k <- 1

for(j in 1:length(psi_vec)){
    for(t in 1:n_reps){
        three_run[[k]] <- restMultiNMTF_run(Xinput = three_data, k_min = 4, 
                                         k_max = 8,
                                psi = (psi_vec[j])*phi_mat, stability = FALSE)
        bisils <- c(sil_score(three_data[[1]], three_run[[k]]$row_clusters[[1]], three_run[[k]]$col_clusters[[1]])$sil, 
                    sil_score(three_data[[2]], three_run[[k]]$row_clusters[[2]], three_run[[k]]$col_clusters[[2]])$sil, 
                    sil_score(three_data[[3]], three_run[[k]]$row_clusters[[3]], three_run[[k]]$col_clusters[[3]])$sil)
        bisils2 <- c(sil_score(three_dt[[1]], three_run[[k]]$col_clusters[[1]], three_run[[k]]$row_clusters[[1]])$sil, 
                    sil_score(three_dt[[2]], three_run[[k]]$col_clusters[[2]], three_run[[k]]$row_clusters[[2]])$sil, 
                    sil_score(three_dt[[3]], three_run[[k]]$col_clusters[[3]], three_run[[k]]$row_clusters[[3]])$sil)
        results[k,] <- c(0, psi_vec[j], jaccard_row(three_run[[k]]$col_clusters[[1]],docs_labs)$f_score[1],
                        jaccard_row(three_run[[k]]$col_clusters[[2]],docs_labs)$f_score[1], 
                        jaccard_row(three_run[[k]]$col_clusters[[3]],docs_labs)$f_score[1], 
                        jaccard_row(three_run[[k]]$col_clusters[[1]],docs_labs)$rec[1],
                        jaccard_row(three_run[[k]]$col_clusters[[2]],docs_labs)$rec[1], 
                        jaccard_row(three_run[[k]]$col_clusters[[3]],docs_labs)$rec[1], 
                        jaccard_row(three_run[[k]]$col_clusters[[1]],docs_labs)$rel[1],
                        jaccard_row(three_run[[k]]$col_clusters[[2]],docs_labs)$rel[1], 
                        jaccard_row(three_run[[k]]$col_clusters[[3]],docs_labs)$rel[1], bisils, bisils2)
         k <- k + 1
  	print(k)
         write.csv(results, "3sources/3sources_psi.csv")
    }
}


res_bbc <- read.csv("3sources/3sources_psi.csv")

#make dataframe 
psi_vec <- seq(0,2000,250)

res_bbc <- res_bbc[,3:(3+4*n_views)]
n_views <- 3
old_names <- c("psi",
                        paste0("F score (V", 1:n_views, ")"),
                 paste0("Recovery (V", 1:n_views,")"),  
                 paste0("Relevance (V", 1:n_views, ")"),
                        paste0("BiS (V", 1:n_views, ")"))
colnames(res_bbc) <- old_names

#calculate means
measures <- c("F score", "BiS", "Recovery", "Relevance")
for(measure in measures){
        m <- paste0(measure,c(" (V1)"," (V2)", " (V3)")) 
        res_bbc <- cbind(res_bbc, (res_bbc[,m[1]]+res_bbc[,m[2]]+res_bbc[,m[3]])/3)
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
write.csv(final_res, "3sources/final_res_resnmtf.csv")
