args <- commandArgs(trailingOnly = TRUE)
dataset <- as.character(args[1])
n_views <- as.numeric(args[2])
# library(ggplot2)
library(latex2exp)
source("SimStudy/Functions/extra_funcs.r")
library(bisilhouette)
if (dataset == "3sources"){
    file_path <- "RealData/3sources/data/three_s_psi_"
    phi_vec <- seq(0, 2000, 50)
}else if (dataset == "single_cell") {
    file_path <- "RealData/single_cell/data/sc_psi_"
    phi_vec <- seq(0, 40000, 1000)
} else if (dataset == "bbcsport") {
    file_path <- "RealData/bbcsport/data/bbc_psi_"
    phi_vec <- seq(0, 2000, 50)
} else {
    file_path <- "RealData/3cancers/data/three_psi_"
    phi_vec <- seq(0, 40000, 1000)
    phi_vec <- c(0, seq(2000, 40000, 1000)) # for 3cancers
}
base_path <- paste0("RealData/", dataset)

n_col <- 3 + 6 * (n_views + 1)

results_euc <- matrix(0, nrow=5*length(phi_vec), ncol=n_col)
results_mann <- matrix(0, nrow=5*length(phi_vec), ncol=n_col)
results_cos <- matrix(0, nrow=5*length(phi_vec), ncol=n_col)

k <- 1
for(phi_val in phi_vec){
    file_cosine <- paste0(file_path, "cosine", phi_val, ".csv")
    file_euc <- paste0(file_path, "euclidean", phi_val, ".csv")
    file_mann <- paste0(file_path, "manhattan", phi_val, ".csv")
    cos_res <- as.matrix(read.csv(file_cosine, row.names=1))
    n_rep_cos <- dim(cos_res)[1] - 1
    euc_res <- as.matrix(read.csv(file_euc, row.names=1))
    n_rep_euc <- dim(euc_res)[1] - 1
    mann_res <- as.matrix(read.csv(file_mann, row.names=1))
    n_rep_mann <- dim(mann_res)[1] - 1
    print(phi_val)
    results_cos[k:(k + n_rep_cos),] <- cos_res
    results_euc[k:(k + n_rep_euc),] <- euc_res
    results_mann[k:(k + n_rep_mann),] <- mann_res
    k <- k + 5
}

write.csv(results_euc, paste0(base_path, "/data/euc_results.csv"))
write.csv(results_cos, paste0(base_path, "/data/cos_results.csv"))
write.csv(results_mann, paste0(base_path, "/data/mann_results.csv"))

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
res_list <- list(results_euc[results_euc[,"psi"]!=0, ], 
            results_cos[results_cos[,"psi"]!=0, ],
            results_mann[results_mann[,"psi"]!=0, ])
dis <- c("E", "C", "M")
sub_res <- vector("list", length=3)
dis_res <- matrix()
for(i in 1:3){
    res <- res_list[[i]]
    # ignore rows with 0 - they didn't converge
    max_bis_e <- which.max(ifelse(res[, "BiS.E"]==0, -Inf, res[, "BiS.E"]))
    max_bis_c <- which.max(ifelse(res[, "BiS.C"]==0, -Inf, res[, "BiS.C"]))
    max_bis_m <- which.max(ifelse(res[, "BiS.M"]==0, -Inf, res[, "BiS.M"]))
    max_f <- which.max(res[, "F.score"])
    sub_res[[i]] <- as.data.frame(cbind(method, res[c(max_bis_e, max_bis_c, max_bis_m, max_f), ]))
    colnames(sub_res[[i]]) <- c("method", old_names)
    write.csv(sub_res[[i]],
            paste0(base_path, "/data/", dis[i], "_results.csv"))
}
dis_study <- rbind(c(sub_res[[1]]$F.score, sub_res[[2]]$F.score, sub_res[[3]]$F.score), 
             c(sub_res[[1]]$psi, sub_res[[2]]$psi, sub_res[[3]]$psi),
             c(sub_res[[1]]$k, sub_res[[2]]$k, sub_res[[3]]$k))
write.csv(dis_study,
    paste0(base_path, "/distance_study.csv"))

#produce plot of f score vs bis 
f_vs_bis_plot <- function(data, dis_name, base_path, col_names, scale=5){
    euc_sc <- read.csv(paste0(base_path, "/data/", dis_name, "_results.csv"), row.names=1)
    colnames(euc_sc) <- col_names

    sc <- as.data.frame(euc_sc)
    sc$BiS <- ifelse(dis_name == "euc", sc$BiS.E, 
                        ifelse(dis_name == "cos", sc$BiS.C, sc$BiS.M))
    p <- ggplot(sc, aes(x = psi)) +
        geom_point(aes(y = F.score, color="F.score"), alpha=0.55) +
        geom_point(aes(y = BiS*scale, color="BiS"), alpha=0.55) + 
        scale_y_continuous(
            name = "F score",
            sec.axis = sec_axis(~ . / 5, name = "BiS") 
        ) +
        labs(x = TeX("$\\phi$")) +
        scale_colour_manual(
            name = "",
            values = c("F.score" = "black", "BiS" = "green"),
            labels = c("F.score" = "F score", "BiS" = "BiS")
        ) +
        theme_minimal() +
        theme(legend.position="none", text = element_text(size = 15))
    
    suppressMessages(ggsave(paste0(base_path, "/f_score_bis_", dis_name, ".pdf"), plot = p, compress = FALSE, device="pdf", width=7,height=7))
}
# f_vs_bis_plot(results_euc, "euc", base_path, old_names)
# f_vs_bis_plot(results_cos, "cos", base_path, old_names, scale=20)
# f_vs_bis_plot(results_mann, "mann", base_path, old_names)


# #bisil plot - fix, need single cell to save row and cols
# if(dataset == "single_cell"){
#     max_f <- which.max(results_euc[, "F.score"])
#     rep_max <- results_euc[max_f,"rep"]
#     psi_max <- results_euc[max_f,"psi"]
#     filepath_row <- paste0(dataset,"/euclidean/data/row_clusts",
#                     psi_max, "_", rep_max, ".xlsx")
#     filepath_col <- paste0(dataset,"/euclidean/data/col_clusts",
#                     psi_max, "_", rep_max, ".xlsx")
#     row_clusts <- import_matrix(filepath_row)
#     col_clusts <- import_matrix(filepath_col)
#     single_cell <- import_matrix("single_cell/data_processed.xlsx")
#     path_to_save <- paste0(base_path,"/sc_bisil_plot.pdf")
#     bisil_plot(single_cell[[1]], row_clusts[[1]], col_clusts[[1]], path_to_save)
# }

