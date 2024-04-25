#bbc overall results
#local paths
filepath <- "3sources"
#source_filepath <- "sim_common_files/simulation_files/"
#source_filepath2 <- "sim_common_files/resNMTF_files/"
#server paths
filepath <- "3sources"
source_filepath <- ""
source_filepath2 <- ""
source(paste0(source_filepath, "gfa_funcs.r"))
source(paste0(source_filepath, "extra_funcs.r"))
source(paste0(source_filepath2, "evaluation_funcs.r"))
source(paste0(source_filepath2, "resNMTF_funcs.r"))
set.seed(40)
three_data <- import_matrix("3sources/3sources_all_diff.xlsx")
docs_labs <- read.csv("3sources/true_labels.csv", row.names=1)

#gfa
n_reps <- 5
n_views <- length(three_data)
bbc_res2 <- vector("list", length=n_reps*2)
n_col <- (n_views+1)*4+2
results <- matrix(0, nrow=n_reps*2, ncol=n_col)
colnames(results) <- c("rep", "method",
                 paste0("F score (V", 1:n_views, ")"), "F score", 
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
                 paste0("BiS-E (V", 1:n_views, ")"), "BiS-E")
k<-1
for(i in 1:n_reps){
        bbc_res2[[k]] <- gfa_apply(lapply(three_data,t), dim(docs_labs)[2])
 
        #assess performance 
        bisils <- calc_all_sils(lapply(three_data,t), bbc_res2[[k]])
        jaccs <-  all_jaccs(docs_labs, bbc_res2[[k]]$row_clusters)
        results[k, ] <- c(i, "gfa", jaccs, bisils$euc)
        k <- k+1
        write.csv(results, paste0(filepath, "/other_method_results.csv"))
}

set.seed(25)
k <- 6
#concatenated data
conc_bbc <- rbind(three_data[[1]], three_data[[2]], three_data[[3]])
n_views <- 1
n_col2 <- (2)*4+2
results <- matrix(0, nrow=n_reps*2, ncol=n_col2)
colnames(results) <- c("rep", "method",
                 paste0("F score (V", 1:n_views, ")"), "F score", 
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
                 paste0("BiS-E (V", 1:n_views, ")"), "BiS-E")
for(t in 1:n_reps){
        bbc_res2[[k]] <- restMultiNMTF_run(Xinput = list(conc_bbc), k_min = 4, 
                                         k_max = 8,stability=FALSE)
        save_results_real(bbc_res2[[k]], paste0(filepath,"/conc_res"),k)
        #assess performance 
        bisils <- calc_all_sils(list(conc_bbc), bbc_res2[[k]])
        jaccs <-  all_jaccs(docs_labs, bbc_res2[[k]]$col_clusters)
        results[k, ] <- c(t, "conc", jaccs, bisils$euc)
        k <- k+1
        write.csv(results, paste0(filepath, "/other_method_results.csv"))
}


results_issvd <- matrix(0, nrow=n_reps, ncol=n_col)
colnames(results_issvd) <- c("rep", "method",
                 paste0("F score (V", 1:n_views, ")"), "F score", 
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
                 paste0("BiS-E (V", 1:n_views, ")"), "BiS-E")
k <- 1
#save results from python doc
for(t in 1:n_reps){
        row_path <- paste0(filepath, "/issvd_res/",t-1,"_row_clusts.xlsx")
        col_path <- paste0(filepath, "/issvd_res/",t-1,"_col_clusts.xlsx")
        row_clusts <- import_matrix(row_path)
        col_clusts <- import_matrix(col_path)
        bisils <- calc_all_sils(three_data, 
                        list("row_clusters" = col_clusts, "col_clusters" = row_clusts))
        jaccs <-  all_jaccs(docs_labs, row_clusts)
        results_issvd[k, ] <- c(t, "issvd", jaccs, bisils$euc)
        k <- k + 1
}
write.csv(results_issvd, paste0(filepath, "/python_method_results.csv"))

# #combine results
filepath <- "3sources"
resnmtf_res <- read.csv(paste0(filepath, "/final_res_resnmtf.csv"))
conc_res <- read.csv(paste0(filepath, "/other_method_results.csv"))[6:10,]
gfa_res <- read.csv(paste0(filepath, "/gfa_method_results.csv"))[1:5,]
issvd_res <- read.csv(paste0(filepath, "/python_method_results.csv"))

colnames <- c("Method",
                 paste0("F score (V", 1:n_views, ")"), "F score", 
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance", 
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery", 
                 paste0("BiS (V", 1:n_views, ")"), "BiS")

gfa_result <- c("GFA",as.matrix(gfa_res[which.max(gfa_res$F.score), 4:ncol(gfa_res)]))
conc_result <- conc_res[which.max(conc_res$F.score),3:ncol(conc_res)]
conc_result <- c("ConcNMTF",rep(conc_result$F.score,n_views+1),
                       rep(conc_result$Relevance,n_views+1),
                       rep(conc_result$Recovery,n_views+1),
                       rep(conc_result$BiS.E,n_views+1))
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