euc <- read.csv("RealData/single_cell/euc_results.csv", row.names = 1)
dataset <- "single_cell"
n_views <- 2
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

print(cor(euc$F.score, euc$BiS.E, method = "spearman"))
print(cor(euc$F.score, euc$Sil, method = "spearman"))

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
three_data <- import_matrix("RealData/3cancers", "tgca_3cancers")
three_data <- lapply(three_data, function(x) x[, 2:1242])
three_data <- lapply(three_data, function(x) apply(x, 2, as.numeric))
docs_labs <- read.csv("RealData/3cancers/tgca_labels.csv")
jacc_mat <- read.csv("RealData/3cancers/data/stability/jacc_mat_5.csv")
jacc_mat < 0.85

rows_1 <- read.csv("RealData/3cancers/data/stability/row_clusts_5_1.csv")
cols_1 <- read.csv("RealData/3cancers/data/stability/col_clusts_5_1.csv")
rows_2 <- read.csv("RealData/3cancers/data/stability/row_clusts_5_2.csv")
cols_2 <- read.csv("RealData/3cancers/data/stability/col_clusts_5_2.csv")
rows_1_sub <- rows_1
rows3_1 <- read.csv("RealData/3cancers/data/stability/row_clusts_3_1.csv")
cols3_1 <- read.csv("RealData/3cancers/data/stability/col_clusts_3_1.csv")
library(bisilhouette)

rows_1_sub <- rows_1
rows_1_sub[, 2:4] <- 0
cols_1_sub <- cols_1
cols_1_sub[, 2:4] <- 0
rows_2_sub <- rows_2
rows_2_sub[, 2:4] <- 0
cols_2_sub <- cols_2
cols_2_sub[, 2:4] <- 0

bisil <- bisilhouette(scale(three_data[[1]]), rows_1_sub, cols_1_sub, "cosine", n_reps = 1)

bisil_2 <- bisilhouette(scale(three_data[[2]]), rows_2_sub, cols_2_sub, "cosine", n_reps = 2)

bisil_3_1 <- bisilhouette(three_data[[1]], rows_1_sub, cols_1_sub, "cosine", n_reps = 2)

bisil_plot(
  scale(three_data[[1]]), rows, cols,
  "RealData/3sources/data/stability/bisilhouette_1_1.pdf", "cosine", 7, 7
)

vals <- bisilhouette(three_data[[1]], rows, cols, "cosine")

cols <- read.csv("RealData/3cancers/data/stability/col_clusts_2_1.csv")
jaccard_row(as.matrix(cols[, c(2, 4)]), as.matrix(docs_labs), print = FALSE)$f_score[1]


jaccard((1:1241)[cols[4] == 1], (1:1241)[docs_labs[3] == 1])
