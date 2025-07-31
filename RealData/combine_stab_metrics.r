args <- commandArgs(trailingOnly = TRUE)
dataset <- as.character(args[1])
n_views <- as.numeric(args[2])
source("SimStudy/Functions/extra_funcs.r")
library(latex2exp)
library(ggplot2)
library(bisilhouette)

file_path <- paste0("RealData/", dataset)
stab_vec <- seq(0, 1, 0.05)
n_col <- 4 + 6 * (n_views + 1)
dataset2 <- ifelse(dataset == "single_cell", "sc", dataset)
dataset2 <- ifelse(dataset == "bbcsport", "bbc", dataset2)
results_bis <- read.csv(
  paste0(file_path, "/", dataset2, "_stab_ResNMTF_BiS.csv"),
  row.names = 1
)
results_fscore <- read.csv(
  paste0(file_path, "/", dataset2, "_stab_ResNMTF_F.csv"),
  row.names = 1
)
results_nmtf <- read.csv(
  paste0(file_path, "/", dataset2, "_stab_NMTF.csv"),
  row.names = 1
)

# process results
old_names <- c(
  "rep", "omega",
  paste0("F score (V", 1:n_views, ")"), "F.score",
  paste0("Relevance (V", 1:n_views, ")"), "Relevance",
  paste0("Recovery (V", 1:n_views, ")"), "Recovery",
  paste0("BiS - E (V", 1:n_views, ")"), "BiS.E",
  paste0("BiS - C (V", 1:n_views, ")"), "BiS.C",
  paste0("BiS - M (V", 1:n_views, ")"), "BiS.M", "k"
)

colnames(results_bis) <- old_names
colnames(results_fscore) <- old_names
colnames(results_nmtf) <- old_names
method <- c(paste0("ResNMTF - ", c("BiS/BiS", "BiS/F", "F/F")), "NMTF")
res_list <- list(
  results_bis,
  results_bis,
  results_fscore,
  results_nmtf
)
sub_res <- vector("list", length = 3)
results <- matrix(0, nrow = 4, ncol = n_col)
for (i in 1:4) {
  res <- res_list[[i]]
  max_bis_e <- which.max(ifelse(res[, "BiS.E"] == 0, -Inf, res[, "BiS.E"]))
  max_bis_c <- which.max(ifelse(res[, "BiS.C"] == 0, -Inf, res[, "BiS.C"]))
  max_bis <- ifelse(
    dataset == "single_cell",
    max_bis_e,
    max_bis_c
  )
  max_f <- which.max(res[, "F.score"])
  row_res <- cbind(
    method[i],
    res[ifelse(i == 1, max_bis, max_f), ][1, ]
  )
  results[i, ] <- unlist(row_res)
}
colnames(results) <- c("method", old_names)
write.csv(
  results,
  paste0(file_path, "/stability_study.csv")
)

if (dataset == "3sources") {
  for (score in c("ResNMTF_BiS", "ResNMTF_F")) {
    df <- read.csv(paste0(file_path, "/", dataset2, "_stab_", score, ".csv"),
      row.names = 1
    )
    colnames(df) <- old_names
    data <- as.data.frame(df) # Create the line plot
    data <- data[data$omega != "none", ]
    data$rep <- factor(data$rep)
    p <- ggplot(data, aes(x = omega, group = rep)) +
      geom_line(aes(y = Relevance), color = "black", alpha = 0.5) +
      scale_x_discrete(breaks = seq(0, 1, 0.1)) +
      # geom_line(aes(y= BiS.E), color="green") +
      labs(
        x = TeX("$\\omega$"), ,
        y = "Relevance"
      ) +
      theme_minimal() +
      theme(text = element_text(size = 15))
    ggsave(paste0(file_path, "/stab_plot_", score, ".pdf"),
      p,
      width = 7, height = 7
    )
  }
  # save bisil plot
  # rep 1 corresponds to the best result in terms of F score
  rows_1 <- read.csv(
    "RealData/3sources/data/stability_ResNMTF_F/row_clusts_1_1.csv"
  )
  cols_1 <- read.csv(
    "RealData/3sources/data/stability_ResNMTF_F/col_clusts_1_1.csv"
  )
  data <- import_matrix("RealData/3sources", "3sources_all_diff")
  bisil_plot(scale(data[[1]]), rows_1, cols_1,
    filename = "RealData/3sources/data/stability_ResNMTF_F/bisil_plot.pdf",
    method = "cosine", w = 7, h = 7
  )
}
