args <- commandArgs(trailingOnly = TRUE)
dis <- as.character(args[1])
# single cell psi selection
source("SimStudy/Functions/evaluation_funcs.r")
source("SimStudy/Functions/extra_funcs.r")
source("visualisation.r")
source("gfa_funcs.r")
if (!requireNamespace("R.matlab", quietly = TRUE)) {
  install.packages("R.matlab")
}
if (!requireNamespace("viridis", quietly = TRUE)) {
  install.packages("viridis")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("latex2exp", quietly = TRUE)) {
  install.packages("latex2exp")
}
library(resnmtf)

labs <- read.csv("RealData/single_cell/true_labs.csv")[, 2:4]
cell_data <- import_matrix("RealData/single_cell/data_processed.xlsx")
phi_mat <- matrix(0, 2, 2)
phi_mat[1, 2] <- 1

n_views <- length(cell_data)
phi_vec <- seq(0, 40000, 500)
n_reps <- 5
n_col <- (n_views + 1) * 6 + 3
results <- matrix(0, nrow = n_reps, ncol = n_col)
colnames(results) <- c(
  "rep", "phi",
  paste0("F score (V", 1:n_views, ")"), "F score",
  paste0("Relevance (V", 1:n_views, ")"), "Relevance",
  paste0("Recovery (V", 1:n_views, ")"), "Recovery",
  paste0("BiS-E (V", 1:n_views, ")"), "BiS-E",
  paste0("BiS-M (V", 1:n_views, ")"), "BiS-M",
  paste0("BiS-C (V", 1:n_views, ")"), "BiS-C", "k"
)
k <- 1

for (phi_val in phi_vec) {
  for (t in 1:n_reps) {
    res <- apply_resnmtf(
      Xinput = cell_data, k_min = 3,
      k_max = 6, phi = phi_val * phi_mat,
      distance = dis, stability = FALSE
    )
    results[k, ] <- c(
      t, phi_val,
      dis_results(cell_data, labs, res, phi_val, t,
        paste0("single_cell/", dis),
        row_same = TRUE
      )
    )
    k <- k + 1
    write.csv(results, paste0("single_cell/data/sc_phi_", dis, phi_val, ".csv"))
  }
}
