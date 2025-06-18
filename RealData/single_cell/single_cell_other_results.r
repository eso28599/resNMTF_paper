filepath <- "RealData/single_cell"
source("SimStudy/OtherMethods/gfa_funcs.r")
source("SimStudy/Functions/evaluation_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
set.seed(40)
labs <- read.csv(paste0(filepath, "/true_labs.csv"))[,2:4]
cell_data <- import_matrix("RealData/single_cell", "data_processed")

#gfa
n_reps <- 5
n_views <- length(cell_data)
bbc_res2 <- vector("list", length=n_reps*2)
n_col <- (n_views+1)*4+2
# results <- matrix(0, nrow=n_reps, ncol=n_col)
# colnames(results) <- c("rep", "method",
#                  paste0("F score (V", 1:n_views, ")"), "F score", 
#                  paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
#                  paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
#                  paste0("BiS-E (V", 1:n_views, ")"), "BiS-E")
# k <- 2
# for(i in 1:n_reps){
#         set.seed(39 + k)
#         bbc_res2[[k]] <- gfa_apply(cell_data, dim(labs)[2])
#         #assess performance 
#         bisils <- calc_all_sils(cell_data, bbc_res2[[k]])
#         jaccs <-  all_jaccs(labs, bbc_res2[[k]]$row_clusters)
#         results[k, ] <- c(i, "gfa", jaccs, bisils$euc)
#         k <- k+1
#         write.csv(results, paste0(filepath, "/gfa_method_results.csv"))
# }

n_views <- length(cell_data)
results_issvd <- matrix(0, nrow=n_reps, ncol=n_col)
colnames(results_issvd) <- c("rep", "method",
                 paste0("F score (V", 1:n_views, ")"), "F score", 
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
                 paste0("BiS-E (V", 1:n_views, ")"), "BiS-E")
k <- 1
#save results from python doc
for(t in 1:n_reps){
        base_path <- paste0(filepath, "/issvd_res")
        row_clusts <- import_matrix(base_path, paste0(t-1,"_row_clusts"))
        col_clusts <- import_matrix(base_path, paste0(t-1,"_col_clusts"))
        bisils <- calc_all_sils(cell_data, 
                        list("row_clusters" = row_clusts, "col_clusters" = col_clusts))
        jaccs <-  all_jaccs(labs, row_clusts)
        results_issvd[k, ] <- c(t, "issvd", jaccs, bisils$euc)
        k <- k + 1
        write.csv(results_issvd, paste0(filepath, "/python_method_results.csv"))
}


n_views <- 2
# #combine results
resnmtf_res <- read.csv(paste0(filepath, "/stability_study.csv"), row.names=1)
gfa_res <- read.csv(paste0(filepath, "/gfa_method_results.csv"))
issvd_res <- read.csv(paste0(filepath, "/python_method_results.csv"))

colnames <- c("Method",
                 paste0("F score (V", 1:n_views, ")"), "F score", 
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
                 paste0("BiS (V", 1:n_views, ")"), "BiS")

gfa_result <- c("GFA",as.matrix(gfa_res[which.max(gfa_res$F.score), 4:ncol(gfa_res)]))
issvd_result <- c("iSSVD", as.matrix(issvd_res[which.max(issvd_res$F.score), 4:ncol(issvd_res)]))
resnmtf_results <- as.matrix(resnmtf_res[,2:ncol(resnmtf_res)])

all_results <- rbind(resnmtf_results, conc_result,
                        gfa_result, issvd_result)
colnames(all_results) <- colnames
all_results <- as.data.frame(all_results)
all_results[,2:ncol(all_results)] <- apply(all_results[,2:ncol(all_results)],2,as.numeric)
write.csv(all_results, paste0(filepath,"/all_results.csv"),row.names=FALSE)
res_text <- kbl(all_results,booktabs=T,"latex",escape=FALSE,digits=4,row.names=FALSE)
df <- kable_styling(res_text)
#save
sink(paste0(filepath, "/all_results.txt"))
print(df)
sink()

colnames_sub <- c("Method", "F score", "Relevance","Recovery", "BiS")
all_results_main <- all_results[,colnames_sub]
write.csv(all_results_main, paste0(filepath,"/all_results_main.csv"),row.names=FALSE)
res_text <- kbl(all_results_main,booktabs=T,"latex",escape=FALSE,digits=4,row.names=FALSE)
df <- kable_styling(res_text)
#save
sink(paste0(filepath, "/all_results_main.txt"))
print(df)
sink()

colnames_sub <- c("Method",  paste0("F score (V", 1:n_views, ")"), "F score")
all_results_main <- all_results[,colnames_sub]
write.csv(all_results_main, paste0(filepath,"/all_results_fscores.csv"),row.names=FALSE)
res_text <- kbl(all_results_main,booktabs=T,"latex",escape=FALSE,digits=4,row.names=FALSE)
df <- kable_styling(res_text)
#save
sink(paste0(filepath, "/all_results_fscores.txt"))
print(df)
sink()
