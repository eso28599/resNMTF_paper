args = commandArgs(trailingOnly = TRUE)
dataset = as.character(args[1])
n_views = as.numeric(args[2])
source("visualisation.r")
source("extra_funcs.r")
library(latex2exp)
library(ggplot2)
if(dataset=="3sources"){
    file_path <- paste0(dataset, "/data/three_s_psi_")
    phi_vec <- c(seq(0, 1800, 50), 1900, 1950)
}else if (dataset == "single_cell") {
    file_path <- paste0(dataset, "/data/sc_phi_")
    phi_vec <- seq(0, 40000, 500)
}else{
    file_path <- paste0(dataset, "/data/", dataset, "_psi_")
    phi_vec <- seq(0, 2000, 50)
}
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

write.csv(results_euc, paste0(dataset, "/data/euc_results.csv"))
write.csv(results_cos, paste0(dataset, "/data/cos_results.csv"))
write.csv(results_mann, paste0(dataset, "/data/mann_results.csv"))

#process results
old_names <- c("rep","psi",
                        paste0("F score (V", 1:n_views, ")"), "F.score",
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance",  
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery",
                        paste0("BiS - E (V", 1:n_views, ")"), "BiS.E",
                        paste0("BiS - C (V", 1:n_views, ")"), "BiS.C", 
                        paste0("BiS - M (V", 1:n_views, ")"), "BiS.M", "k")

colnames(results_cos) <- old_names
colnames(results_euc) <- old_names
colnames(results_mann) <- old_names
method <- paste0("ResNMTF - ", c("BiS (E)", "BiS (C)", "BiS (M)", "F score"))
res_list <- list(results_euc, results_cos, results_mann)
dis <- c("E", "C", "M")
sub_res <- vector("list", length=3)
dis_res <- matrix()
for(i in 1:3){
    res <- res_list[[i]]
    max_bis_e <- which.max(res[, "BiS.E"])
    max_bis_c <- which.max(res[, "BiS.C"])
    max_bis_m <- which.max(res[, "BiS.M"])
    max_f <- which.max(res[, "F.score"])
    sub_res[[i]] <- as.data.frame(cbind(method, res[c(max_bis_e, max_bis_c, max_bis_m, max_f), ]))
    colnames(sub_res[[i]]) <- c("method", old_names)
    write.csv(sub_res[[i]],
            paste0(dataset, "/data/", dis[i], "_results.csv"))
}
dis_study <- rbind(c(sub_res[[1]]$F.score, sub_res[[2]]$F.score, sub_res[[3]]$F.score), 
             c(sub_res[[1]]$psi, sub_res[[2]]$psi, sub_res[[3]]$psi))
write.csv(dis_study,
    paste0(dataset, "/distance_study.csv"))

#produce plot of f score vs bis 
euc_sc <- read.csv(paste0(dataset, "/data/euc_results.csv"), row.names=1)
colnames(euc_sc) <- old_names
sc <- as.data.frame(euc_sc)

p <- ggplot(sc, aes(x = psi)) +
  geom_point(aes(y = F.score, color="F.score"),alpha=0.55) +
  geom_point(aes(y = BiS.E*5,color="BiS.E"),alpha=0.55) + 
  scale_y_continuous(
    name = "F score",
    sec.axis = sec_axis(~ . / 5, name = "BiS") 
  ) +
  labs(x = TeX("$\\phi$")) +
  scale_colour_manual(
    name = "",
    values = c("F.score" = "black", "BiS.E" = "green"),
    labels = c("F.score" = "F score", "BiS.E" = "BiS")
  ) +
  theme_minimal()

suppressMessages(ggsave(paste0(dataset, "_f_score_bis_psi.pdf"), plot = p, compress = FALSE, device="pdf"))


#bisil plot
if(dataset=="single_cell"){
    max_f <- which.max(results_euc[, "F.score"])
    rep_max <- results_euc[max_f,"rep"]
    psi_max <- results_euc[max_f,"psi"]
    filepath_row <- paste0(dataset,"/euclidean/data/row_clusts",
                    psi_max, "_", rep_max, ".xlsx")
    filepath_col <- paste0(dataset,"/euclidean/data/col_clusts",
                    psi_max, "_", rep_max, ".xlsx")
    row_clusts <- import_matrix(filepath_row)
    col_clusts <- import_matrix(filepath_col)
    single_cell <- import_matrix("single_cell/data_processed.xlsx")
    path_to_save <- paste0(dataset,"/sc_bisil_plot.pdf")
    bisil_plot(single_cell[[1]], row_clusts[[1]], col_clusts[[1]], path_to_save)
}

