euc <- read.csv("RealData/single_cell/euc_results.csv", row.names = 1)
dataset <- "single_cell" # or "text"
n_views <- 2 # Number of views, adjust as necessar
old_names <- c(
  "rep", "psi",
  paste0("F score (V", 1:n_views, ")"), "F.score",
  paste0("Relevance (V", 1:n_views, ")"), "Relevance",
  paste0("Recovery (V", 1:n_views, ")"), "Recovery",
  paste0("BiS - E (V", 1:n_views, ")"), "BiS.E",
  paste0("BiS - C (V", 1:n_views, ")"), "BiS.C",
  paste0("BiS - M (V", 1:n_views, ")"), "BiS.M", "k"
)
if (dataset == "single_cell") {
  col_names_euc <- c(
    old_names,
    paste0("Sil (V", 1:n_views, ")"), "Sil"
  )
} else {
  col_names_euc <- old_names
}

cor(euc$F.score, euc$BiS.E, method = "spearman")
cor(euc$F.score, euc$Sil, method = "spearman")

colnames(euc) <- col_names_euc

sc <- as.data.frame(euc)

p <- ggplot(sc, aes(x = psi)) +
  geom_point(aes(y = F.score, color = factor("F.score", levels = c("F.score", "BiS", "Sil"))), alpha = 0.55) +
  geom_point(aes(y = 10 * (BiS.E + 0.06), color = "BiS"), alpha = 0.55) +
  geom_point(aes(y = 10 * (Sil + 0.06), color = "Sil"), alpha = 0.55) +
  scale_y_continuous(
    name = "",
    limits = c(0, 0.6),
    sec.axis = sec_axis(~ . / 10 - 0.06, name = "")
  ) +
  labs(x = TeX("$\\phi$")) +
  scale_colour_manual(
    name = "",
    values = c("F.score" = "black", "BiS" = "green", "Sil" = "blue"),
    labels = c(
      "F.score" = " F score (left axis)", "BiS" = "Bisilhouette (right axis)", "Sil" = "Silhouette (right axis)",
      breaks = c("F.score", "BiS", "Sil")
    )
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", text = element_text(size = 15))

suppressMessages(ggsave(paste0(base_path, "/f_score_bis_sil.pdf"), plot = p, compress = FALSE, device = "pdf", width = 8, height = 8))

source("SimStudy/Functions/extra_funcs.r")
source("SimStudy/Functions/evaluation_funcs.r")
three_data <- import_matrix("RealData/3sources", "3sources_all_diff")
docs_labs <- read.csv("RealData/3cancers/tgca_labels.csv")
rows <- read.csv("RealData/3cancers/data/stability/row_clusts_1_1.csv")
library(bisilhouette)

bisil_plot(
  scale(three_data[[1]]), rows, cols,
  "RealData/3sources/data/stability/bisilhouette_1_1.pdf", "cosine", 7, 7
)

vals <- bisilhouette(three_data[[1]], rows, cols, "cosine")

cols <- read.csv("RealData/3cancers/data/stability/col_clusts_2_1.csv")
jaccard_row(as.matrix(cols[, c(2, 4)]), as.matrix(docs_labs), print = FALSE)$f_score[1]


jaccard((1:1241)[cols[4] == 1], (1:1241)[docs_labs[3] == 1])
