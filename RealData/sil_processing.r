args <- commandArgs(trailingOnly = TRUE)
dataset <- as.character(args[1])
n_views <- as.numeric(args[2])
source("extra_funcs.r")
library(latex2exp)
library(ggplot2)

file_path <- paste0(dataset, "/data/sil_psi_")
phi_vec <- seq(0, 40000, 500)

n_col <- 3 + 6 * (n_views + 1)

results_euc <- matrix(0, nrow = 5 * length(phi_vec), ncol = n_col)

k <- 1
for (phi_val in phi_vec) {
  file_euc <- paste0(file_path, "euclidean", phi_val, ".csv")
  results_euc[k:(k + 4), ] <- as.matrix(read.csv(file_euc, row.names = 1))
  k <- k + 5
}

write.csv(results_euc, paste0(dataset, "/data/sil_results.csv"))

# process results
old_names <- c(
  "rep", "psi",
  paste0("F score (V", 1:n_views, ")"), "F.score",
  paste0("Relevance (V", 1:n_views, ")"), "Relevance",
  paste0("Recovery (V", 1:n_views, ")"), "Recovery",
  paste0("BiS - E (V", 1:n_views, ")"), "BiS.E",
  paste0("BiS - C (V", 1:n_views, ")"), "BiS.C",
  paste0("BiS - M (V", 1:n_views, ")"), "BiS.M", "k"
)

colnames(results_euc) <- old_names
sc <- as.data.frame(results_euc)

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
  theme(legend.position = "none")

suppressMessages(ggsave(paste0(dataset, "_f_score_sil_psi.pdf"),
  plot = p, compress = FALSE, device = "pdf"
))
