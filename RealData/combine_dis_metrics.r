args <- commandArgs(trailingOnly = TRUE)
dataset <- as.character(args[1])
n_views <- as.numeric(args[2])
library(ggplot2)
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
n_col_euc <- ifelse(dataset == "single_cell", n_col + 3, n_col)

results_euc <- matrix(0, nrow=5*length(phi_vec), ncol=n_col_euc)
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
if (dataset == "single_cell") {
    col_names_euc <- c(old_names,
                        paste0("Sil (V", 1:n_views, ")"), "Sil")
} else {
    col_names_euc <- old_names
}

colnames(results_cos) <- old_names
colnames(results_euc) <- col_names_euc
colnames(results_mann) <- old_names
method <- paste0("ResNMTF - ", c("BiS (E)", "BiS (C)", "BiS (M)", "F score"))
res_list <- list(results_euc[results_euc[,"psi"]!=0, 1:n_col], 
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

if (dataset == "single_cell") {
    euc_sc <- read.csv(paste0(base_path, "/data/euc_results.csv"), row.names=1)
    colnames(euc_sc) <- col_names_euc

    sc <- as.data.frame(euc_sc)
    p <- ggplot(sc, aes(x = psi)) +
        geom_point(aes(y = F.score, color=factor("F.score", levels = c("F.score", "BiS", "Sil"))), alpha=0.55) +
        geom_point(aes(y = 10 * (BiS.E + 0.06), color="BiS"), alpha=0.55) + 
        geom_point(aes(y = 10 * (Sil + 0.06), color="Sil"), alpha=0.55) +
        scale_y_continuous(
            name = "",
            limits = c(0, 0.6),
            sec.axis = sec_axis(~ . / 10 - 0.06, name = "") 
        ) +
        labs(x = TeX("$\\phi$")) +
        scale_colour_manual(
            name = "",
            values = c("F.score" = "black", "BiS" = "green", "Sil" = "blue"),
            labels = c("F.score" = " F score (left axis)", "BiS" = "Bisilhouette (right axis)", "Sil" = "Silhouette (right axis)", 
            breaks = c("F.score", "BiS", "Sil") )
        ) +
        theme_minimal() +
        theme(legend.position="bottom", text = element_text(size = 15))
    
    suppressMessages(ggsave(paste0(base_path, "/f_score_bis_sil.pdf"), plot = p, compress = FALSE, device="pdf", width=7,height=8))
}