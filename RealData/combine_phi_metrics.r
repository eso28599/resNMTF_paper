args <- commandArgs(trailingOnly = TRUE)
dataset <- as.character(args[1])
n_views <- as.numeric(args[2])
source("../SimStudy/Functions/extra_funcs.r")
library(latex2exp)
library(ggplot2)
library(bisilhouette)
if (dataset == "3sources") {
  file_path <- paste0(dataset, "/data/three_s_psi_")
  # phi_vec <- c(seq(0, 1800, 50), 1900, 1950) # full distance study
  phi_vec <- c(350, 1700) # end results
} else if (dataset == "single_cell") {
  file_path <- paste0(dataset, "/data/sc_phi_")
  phi_vec <- seq(0, 40000, 500)
} else {
  file_path <- paste0(dataset, "/data/", dataset, "_psi_")
  phi_vec <- seq(0, 2000, 50)
  phi_vec <- c(1850, 600) # end results
}

n_col <- 3 + 4 * (n_views + 1)

results_euc <- matrix(0, nrow = 5 * length(phi_vec), ncol = n_col)

k <- 1
for (phi_val in phi_vec) {
  file_euc <- paste0(file_path, "euclidean", phi_val, ".csv")
  results_euc[k:(k + 4), ] <- as.matrix(read.csv(file_euc, row.names = 1))
  k <- k + 5
}
write.csv(results_euc, paste0(dataset, "/data/euc_results.csv"))

# process results
old_names <- c(
  "rep", "psi",
  paste0("F score (V", 1:n_views, ")"), "F.score",
  paste0("Relevance (V", 1:n_views, ")"), "Relevance",
  paste0("Recovery (V", 1:n_views, ")"), "Recovery",
  paste0("BiS - E (V", 1:n_views, ")"), "BiS.E",
  "k"
)
colnames(results_euc) <- old_names

# produce plot of f score vs bis
euc_sc <- read.csv(paste0(dataset, "/data/euc_results.csv"), row.names = 1)
colnames(euc_sc) <- old_names
sc <- as.data.frame(euc_sc)

p <- ggplot(sc, aes(x = psi)) +
  geom_point(aes(y = F.score, color = "F.score"), alpha = 0.55) +
  geom_point(aes(y = BiS.E * 5, color = "BiS.E"), alpha = 0.55) +
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
  theme_minimal() +
  theme(legend.position = "none", text = element_text(size = 15))

suppressMessages(ggsave(paste0(dataset, "_f_score_bis_psi.pdf"),
  plot = p, compress = FALSE, device = "pdf", width = 7, height = 7
))


# bisil plot
if (dataset == "single_cell") {
  max_f <- which.max(results_euc[, "F.score"])
  rep_max <- results_euc[max_f, "rep"]
  psi_max <- results_euc[max_f, "psi"]
  filepath_row <- paste0(
    dataset, "/euclidean/data/row_clusts",
    psi_max, "_", rep_max, ".xlsx"
  )
  filepath_col <- paste0(
    dataset, "/euclidean/data/col_clusts",
    psi_max, "_", rep_max, ".xlsx"
  )
  row_clusts <- import_matrix(filepath_row)
  col_clusts <- import_matrix(filepath_col)
  single_cell <- import_matrix("single_cell/data_processed.xlsx")
  path_to_save <- paste0(dataset, "/sc_bisil_plot.pdf")
  bisil_plot(single_cell[[1]], row_clusts[[1]], col_clusts[[1]], path_to_save)
}
