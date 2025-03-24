source("SimStudy/RunSim/Functions/evaluation_funcs.r")
library(ggplot2)
library(ggfortify)
library(ggpubr)

break_func <- function(lower) {
  if (min(lower) < 0) {
    return(seq(-1, 1, 0.2)[seq(-1, 1, 0.2) > (min(lower) - 0.2)])
  } else {
    return(seq(0, 1, 0.2))
  }
}

df_plot <- function(values){
  x <- c()
  y <- c()
  cluster <- c()
  for (i in seq_along(length(values))) {
    x <- c(x, seq_along(length(values[[i]])))
    y <- c(y, sort(values[[i]]))
    cluster <- c(cluster, rep(i, length(values[[i]])))
  }
  x <- seq_along(length(x))
  return(data.frame(x = x, y = y, clust = cluster))
}

x_breaks <- function(df) {  
  n_c <- length(unique(df$clust))
  ticks <- c()
  num <- c(0)
  for (i in 1:n_c) {
    num <- c(num, sum(df$clust == i))
    ticks <- c(ticks, floor(sum(df$clust == i) / 2))
  }
  return(cumsum(num[1:(n_c)]) + ticks)
}

bisil_plot <- function(data, rows, cols, filename = NULL) {
    scores <- sil_score_inner(data, rows, cols)
    df <- df_plot(scores$vals)
    break_vals <- break_func(df$y)
    breaks_x <- x_breaks(df)
    bis_val <- scores$sil
    p <- ggplot(df) +
        geom_col(aes(x = x, y = y, group = clust, color = clust, fill = clust), width = 1, position = position_dodge()) +
        scale_y_continuous(breaks = break_vals, limits = c(min(break_vals), max(break_vals))) +
        scale_x_continuous(breaks = breaks_x, labels = 1:(dim(rows)[2])) +
        scale_color_viridis_c() +
        scale_fill_viridis_c() +
        ylab("Bisilhouette score") +
        xlab("Bicluster") +
        geom_hline(yintercept = bis_val, linetype = "dashed", color = "black") +
        coord_flip() +
        theme_minimal() +
        theme(
            legend.position = "none", axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            axis.ticks = element_line(size = 0.4),
            axis.line = element_line(size = 0.4)
        )

    if (!is.null(filename)) {
        ggsave(filename, p)
    }
    return(p)
}

pca_plots <- function(data, rows, cols, filename = NULL) {
    n_c <- dim(rows[[1]])[2]
    n_v <- length(data)
    plots <- vector("list", length = n_c * n_v)
    for (view in 1:n_v) {
        clust_vector <- apply(rows[[view]], 1, function(x) which(x == 1, arr.ind = TRUE))
        for (i in 1:n_c) {
            clust_data <- data[[view]][, cols[[view]][, i] == 1]
            pca_res <- prcomp(clust_data, scale. = TRUE)
            color_vec <- unlist(lapply(clust_vector, function(x) ifelse(length(x) == 0, 10, ifelse(any(x == i), i, x[1]))))
            all_data <- data.frame(clust_data, "Bicluster" = color_vec)
            all_data$Bicluster <- factor(all_data$Bicluster)
            p <- autoplot(pca_res, data = all_data, colour = "Bicluster", legend = "bottom") +
                scale_colour_viridis_d() +
                theme_minimal(base_size = 9) +
                theme(
                    legend.position = "bottom",
                    axis.ticks = element_line(size = 0.3),
                    axis.line = element_line(size = 0.3)
                )
            plots[[n_v * (i - 1) + view]] <- p
        }
    }
    p <- ggarrange(
        plotlist = plots,
        labels = "AUTO",
        ncol = n_v, nrow = n_c, common.legend = TRUE, legend = "top", font.label = list(size = 10)
    )
    if (!is.null(filename)) {
        ggsave(filename, p)
    }
    return(p)
}
