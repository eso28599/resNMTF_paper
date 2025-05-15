source("SimStudy/Functions/data_generation.r")
source("SimStudy/Functions/extra_funcs.r")
# load libraries
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(resnmtf)

phi <- matrix(0, nrow = 2, ncol = 2)
phi[1, c(2)] <- 200
set.seed(15)
save_data(
  rep(200, 2), c(100, 250),
  4,
  "Exploration", 3, TRUE, FALSE
)

data4 <- import_matrix("Exploration/data.xlsx")
nmtf_results4 <- apply_resnmtf(Xinput = data4, phi = phi)

# shuffle data
updated_data <- vector("list", length = 2)
updated_data[[1]] <- matrix(sample(data4[[1]]), 200, 100)
updated_data[[2]] <- matrix(sample(data4[[2]]), 200, 250)
nmtf_shuffled <- apply_resnmtf(
  Xinput = updated_data,
  phi = phi, KK = rep(4, 2)
)
# plot distributions
n_views <- 2
n_bicls <- 4
plots_true <- vector("list", length = n_bicls)
plots_shuff <- vector("list", length = n_bicls)
n <- 1
for (i in 1:1) {
  for (j in 1:n_bicls) {
    dens1 <- density(nmtf_results4$Foutput[[i]][, j],
      from = 0, to = max(nmtf_results4$Foutput[[i]])
    )
    dens2 <- density(nmtf_shuffled$Foutput[[i]][, j],
      from = 0, to = max(nmtf_results4$Foutput[[i]])
    )
    m_fac1 <- 1 +
      sum(dens1$y[dens1$x > max(nmtf_results4$Foutput[[i]][, j])]) /
        sum(dens1$y[dens1$x < max(nmtf_results4$Foutput[[i]][, j])])
    m_fac2 <- 1 +
      sum(dens2$y[dens2$x > max(nmtf_shuffled$Foutput[[i]][, j])]) /
        sum(dens2$y[dens2$x < max(nmtf_shuffled$Foutput[[i]][, j])])
    dens1$y[dens1$x > max(nmtf_results4$Foutput[[i]][, j])] <- 0
    dens2$y[dens2$x > max(nmtf_shuffled$Foutput[[i]][, j])] <- 0
    dens1$y <- (dens1$y) * m_fac1
    dens2$y <- (dens2$y) * m_fac2
    df1 <- data.frame("Entry" = dens1$x, "Density" = dens1$y)
    df2 <- data.frame("Entry" = dens2$x, "Density" = dens2$y)
    xlims <- c(0, 0.032)
    ylims <- c(0, 360)

    if ((j == 1) | (j == 2)) {
      plots_true[[n]] <- ggplot(df1, aes(x = Entry, y = Density)) +
        geom_line() +
        xlim(xlims[1], xlims[2]) +
        ylim(ylims[1], ylims[2]) +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(), # remove x axis labels
          axis.title.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth = 0.3),
          axis.line = element_line(size = 0.3)
        )
      plots_shuff[[n]] <- ggplot(df2, aes(x = Entry, y = Density)) +
        geom_line() +
        xlim(xlims[1], xlims[2]) +
        ylim(ylims[1], ylims[2]) +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(), # remove x axis labels
          axis.title.y = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(linewidth = 0.3),
          axis.line = element_line(size = 0.3)
        )
    } else {
      plots_true[[n]] <- ggplot(df1, aes(x = Entry, y = Density)) +
        geom_line() +
        xlim(xlims[1], xlims[2]) +
        ylim(ylims[1], ylims[2]) +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(), # remove x axis labels
          axis.title.y = element_blank(),
          axis.ticks.y = element_line(linewidth = 0.3),
          axis.line = element_line(linewidth = 0.3)
        )
      plots_shuff[[n]] <- ggplot(df2, aes(x = Entry, y = Density)) +
        geom_line() +
        xlim(xlims[1], xlims[2]) +
        ylim(ylims[1], ylims[2]) +
        theme_minimal() +
        theme(
          axis.title.x = element_blank(), # remove x axis labels
          axis.title.y = element_blank(),
          axis.ticks.y = element_line(linewdith = 0.3),
          axis.line = element_line(size = 0.3)
        )
    }

    n <- n + 1
  }
}

plots_true[[2]] <- plots_true[[2]] + rremove("y.text") + rremove("y.ticks")
plots_true[[4]] <- plots_true[[4]] + rremove("y.text") + rremove("y.ticks")
plots_shuff[[2]] <- plots_shuff[[2]] + rremove("y.text") + rremove("y.ticks")
plots_shuff[[4]] <- plots_shuff[[4]] + rremove("y.text") + rremove("y.ticks")
figure <- ggarrange(
  plotlist = plots_true[1:4], # remove axis labels from plots
  labels = NULL,
  ncol = 2, nrow = 2,
  align = "hv",
  font.label = list(
    size = 10, color = "black", face = "bold",
    family = NULL, position = "top"
  )
)

figure <- annotate_figure(figure,
  left = textGrob("Density", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
  bottom = textGrob("F entry", gp = gpar(cex = 1.3))
)

figure_shuff <- ggarrange(
  plotlist = plots_shuff[1:4], # remove axis labels from plots
  labels = NULL,
  ncol = 2, nrow = 2,
  align = "hv",
  font.label = list(
    size = 10, color = "black",
    face = "bold", family = NULL, position = "top"
  )
)

figure_shuff <- annotate_figure(figure_shuff,
  left = textGrob("Density", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
  bottom = textGrob("F entry", gp = gpar(cex = 1.3))
)
suppressMessages(ggsave("Exploration/f_data_dists.pdf",
  plot = figure, compress = FALSE,
  device = "pdf", width = 18, height = 22, units = "cm"
))
suppressMessages(ggsave("Exploration/shuffled_data_dists.pdf",
  plot = figure_shuff, compress = FALSE,
  device = "pdf", width = 18, height = 22, units = "cm"
))
