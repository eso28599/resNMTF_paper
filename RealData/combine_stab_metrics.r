args = commandArgs(trailingOnly = TRUE)
dataset = as.character(args[1])
n_views = as.numeric(args[2])
source("visualisation.r")
source("extra_funcs.r")
library(latex2exp)
library(ggplot2)

file_path <- paste0(dataset, "/data/", dataset, "_psi_")
stab_vec <- seq(0, 1, 0.05)

n_col <- 4 + 6 * (n_views + 1)

dataset2 <- ifelse(dataset=="single_cell", "sc", dataset)
results_bis <- read.csv(paste0(dataset, "/", dataset2, "_stab_bis.csv"), row.names=1)
results_fscore <- read.csv(paste0(dataset, "/", dataset2, "_stab_fscore.csv"), row.names=1)
results_nmtf <- read.csv(paste0(dataset, "/", dataset2, "_stab_nmtf.csv"), row.names=1)


#process results
old_names <- c("rep","omega",
                        paste0("F score (V", 1:n_views, ")"), "F.score",
                 paste0("Relevance (V", 1:n_views, ")"), "Relevance",  
                 paste0("Recovery (V", 1:n_views, ")"), "Recovery",
                        paste0("BiS - E (V", 1:n_views, ")"), "BiS.E",
                        paste0("BiS - C (V", 1:n_views, ")"), "BiS.C", 
                        paste0("BiS - M (V", 1:n_views, ")"), "BiS.M", "k")

colnames(results_bis) <- old_names
colnames(results_fscore) <- old_names
colnames(results_nmtf) <- old_names
method <- paste0("ResNMTF - ", c("BiS (S)", "F score (S)"))
res_list <- list(results_bis,
            results_fscore,
            results_nmtf)
sub_res <- vector("list", length=3)
results <- matrix(0, nrow=3, ncol=n_col)
for(i in 1:3){
    res <- res_list[[i]]
    max_bis_e <- which.max(res[, "BiS.E"])
    max_f <- which.max(res[, "F.score"])
    sub_res[[i]] <- as.data.frame(cbind(method[i], res[ifelse(i==1, max_bis_e, max_f), ]))
    colnames(sub_res[[i]]) <- c("method", old_names)
}
dis_study <- rbind(c(sub_res[[1]]$F.score, sub_res[[2]]$F.score, sub_res[[3]]$F.score), 
             c(sub_res[[1]]$omega, sub_res[[2]]$omega, sub_res[[3]]$omega),
             c(sub_res[[1]]$k, sub_res[[2]]$k, sub_res[[3]]$k))
write.csv(dis_study,
    paste0(dataset, "/stability_study.csv"))

if(dataset=="bbc"){
    for(score in c("bis", "fscore")){
        df <- read.csv(paste0(dataset, "/", dataset2, "_stab_",score,".csv"), row.names=1)
        data <- as.data.frame(df)# Create the line plot
        data <- data[data$omega!="none",]
        data$rep <- factor(data$rep)
        p <- ggplot(data, aes(x = omega, group = rep)) +
        geom_line(aes(y= F.score), color="black", alpha=0.5) +
        geom_line(aes(y= BiS.E), color="green") +
        labs(
            x = TeX("$\\omega$"),,
            y = "Measure") +
        theme_minimal()+
        theme(text = element_text(size = 11))
        ggsave(paste0(dataset,"/stab_plot_", score, ".pdf"),p, width=7,height=7)
    }
    
}
