# extracting scores
k_val <- 10
results <- apply_resnmtf(data,
  phi = phi, k_val = k_val,
  stability = FALSE,
  spurious = FALSE
)

library(ggplot2)
library(reshape2)
# data <- # your list of data matrices
n_views <- 4
view_names <- c("mrna", "protein", "chromatin", "something") # idk
k_val <- 3 # number of factors to consider
x_axis_label <- "View" # label for x-axis - could be "omic", or blank ("")
scaled_data <- lapply(data, function(x) apply(x, 2, function(x) x / sum(x)))
factor <- "both"
r2_scores <- matrix(0, nrow = n_views, ncol = k_val)
r2_scores_per_view <- numeric(n_views)
for (n in 1:n_views) {
  b <- sum(scaled_data[[n]]**2)
  for (k in 1:k_val) {
    # row factor k alone
    if (factor == "rows") {
      recon <- matrix(results$output_f[[n]][, k]) %*%
        ((results$output_s[[n]] %*% t(results$output_g[[n]]))[k, ])
    } else if (factor == "cols") {
      recon <- (results$output_f[[n]] %*% results$output_s[[n]])[, k] %*%
        t(results$output_g[[n]][, k])
    } else {
      recon <- results$output_s[[n]][k, k] * results$output_f[[n]][, k] %*%
        t(results$output_g[[n]][, k])
    }
    a <- sum((scaled_data[[n]] - recon)**2)
    r2_scores[n, k] <- 1 - a / b
  }
  recon_vn <- (results$output_f[[n]] %*% results$output_s[[n]]) %*%
    t(results$output_g[[n]])
  a <- sum((scaled_data[[n]] - recon_vn)**2)
  r2_scores_per_view[n] <- 1 - a / b
}
rownames(r2_scores) <- view_names
colnames(r2_scores) <- 1:k_val
r2_mk_df <- melt(
  r2_scores * 100,
  varnames = c("View", "Factor"), value.name = "value"
)
r2_mk_df$View <- factor(r2_mk_df$View, levels = view_names)
r2_mk_df$Factor <- factor(r2_mk_df$Factor, levels = 1:k_val)
r2_m_df <- data.frame(
  View = factor(view_names, levels = view_names),
  R2 = r2_scores_per_view * 100
)

p1 <- ggplot(r2_mk_df, aes(x = View, y = Factor)) +
  geom_tile(aes(fill = value), color = "black") +
  scale_fill_viridis_c(
    guide = "colorbar",
    limits = c(min(r2_scores * 100), max(r2_scores * 100))
  ) +
  labs(x = x_axis_label, y = "Factor", title = "") +
  guides(fill = guide_colorbar("Var. (%)")) +
  theme(
    axis.text.x = element_text(size = rel(1.0), color = "black"),
    axis.text.y = element_text(size = rel(1.1), color = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = rel(1.0))
  )
viridis_color <- "#440154FF"
p2 <- ggplot(r2_m_df, aes(x = View, y = R2)) +
  geom_bar(
    stat = "identity",
    fill = viridis_color, colour = viridis_color, width = 0.9
  ) +
  xlab(x_axis_label) +
  scale_fill_viridis_c(guide = "colorbar") +
  ylab("Variance explained (%)") +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    panel.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text()
  )

results$output_f[[1]]
