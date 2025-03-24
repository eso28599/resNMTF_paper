suppressWarnings(suppressPackageStartupMessages(library("rio")))
library(ggplot2)
library("tidyr")
library("dplyr")
library(latex2exp)
library(viridis)
# import as list of matrices instead of as a list
import_matrix <- function(filename) {
  return(lapply(import_list(filename), function(x) as.matrix(x)))
}
calc_averages <- function(mat_res) {
  return(c(rowMeans(mat_res), mean(mat_res)))
}

get_averages <- function(list_avgs, data_name, filename, repeat_no, measure) {
  results <- import_matrix(paste0(
    paste0(data_name, "/"),
    paste0(filename, ".xlsx")
  ))
  for (j in seq_len(length(measure))) {
    list_avgs[[j]][repeat_no, ] <- mean(results[[measure[j]]])
  }
  return(list_avgs)
}


avg_over_repeats <- function(avgs, i, measure, name1, name2, cluster) {
  rep_avgs <- colMeans(avgs[[i]])
  rep_sd <- apply(avgs[[i]], 2, sd)
  return(cbind(name1, name2, cluster, measure, rep_avgs, rep_sd))
}

make_plot_line <- function(data_frame, cluster, plot_title, sub_title,
                           x_title, y_title, filename, measures) {
  data_frame[["Measure"]] <- factor(data_frame[["Measure"]],
    labels = c("F score", "Relevance", "Recovery", "CSR", "BiS")
  )
  if (min(data_frame$Mean) < 0) {
    y_min <- min(data_frame$Mean)
    y_lims <- c(y_min, 1)
    if (y_min < -0.2) {
      y_breaks <- c(
        sort(seq(-0.2, y_min, by = -0.2)),
        seq(0, 1, length = 6)
      )
    } else {
      y_breaks <- c(-0.2, seq(0, 1, length = 6))
    }
  } else {
    y_lims <- c(0, 1)
    y_breaks <- seq(0, 1, length = 6)
  }
  if (x_title == "Sigma") {
    x_title <- TeX("$\\sigma$")
  }
  x_breaks <- unique(data_frame$Factor)
  p <- ggplot(subset(data_frame, Choice == cluster), aes(
    x = Factor, y = Mean,
    group = Measure, color = Measure, linetype = Measure
  )) +
    geom_line(stat = "identity") +
    geom_point(show.legend = FALSE) +
    geom_errorbar(aes(ymin = Mean - S.d., ymax = Mean + S.d.),
      width = .2,
      position = position_dodge(0.01), show.legend = FALSE, linetype = 1
    ) +
    expand_limits(y = y_lims) +
    scale_color_viridis_d() +
    scale_y_continuous(breaks = y_breaks) +
    scale_x_discrete(breaks = x_breaks) +
    labs(
      title = plot_title,
      subtitle = sub_title,
      x = x_title,
      y = y_title
    )
  p <- p +
    theme_minimal() +
    theme(
      text = element_text(size = 11), legend.position = "bottom",
      axis.ticks = element_line(size = 0.3),
      axis.line = element_line(size = 0.3)
    )
  suppressMessages(ggsave(filename, plot = p, compress = FALSE, device = "pdf"))
  return(p)
}

make_plot_measure <- function(data_frame, cluster, plot_title, sub_title,
                              x_title, y_title, filename, measure, methods) {
  data_frame <- subset(data_frame, Measure == measure)
  if (min(data_frame$Mean) < 0) {
    y_min <- min(data_frame$Mean)
    y_lims <- c(y_min, 1)
    if (y_min < -0.2) {
      y_breaks <- c(
        sort(seq(-0.2, y_min, by = -0.2)),
        seq(0, 1, length = 6)
      )
    } else {
      y_breaks <- c(-0.2, seq(0, 1, length = 6))
    }
  } else {
    y_lims <- c(0, 1)
    y_breaks <- seq(0, 1, length = 6)
  }
  if (x_title == "Sigma") {
    x_title <- TeX("$\\sigma$")
  }
  p <- ggplot(subset(data_frame, Choice == cluster), aes(
    x = Factor, y = Mean,
    group = Method, color = Method, linetype = Method
  )) +
    geom_line(stat = "identity") +
    geom_point(show.legend = FALSE) +
    geom_errorbar(aes(ymin = Mean - S.d., ymax = Mean + S.d.),
      width = .2,
      position = position_dodge(0.01), linetype = 1,
      show.legend = FALSE
    ) +
    expand_limits(y = y_lims) +
    scale_color_viridis_d() +
    scale_y_continuous(breaks = y_breaks)
  p <- p +
    theme_minimal() +
    theme(
      text = element_text(size = 11), legend.position = "none",
      axis.ticks = element_line(size = 0.3),
      axis.line = element_line(size = 0.3)
    )
  p_save <- p + labs(x = x_title, y = y_title)
  p <- p + labs(
    x = x_title,
    y = y_title
  )

  suppressMessages(ggsave(filename, plot = p, compress = FALSE, device = "pdf"))
  return(list("a" = p_save, "b" = p))
}

make_plot_line_phi <- function(data_frame, cluster, plot_title,
                               sub_title, x_title, y_title,
                               filename, measures) {
  data_frame[["Measure"]] <- factor(data_frame[["Measure"]],
    labels = measures
  )

  if (min(data_frame$Mean) < 0) {
    y_min <- min(data_frame$Mean)
    y_lims <- c(y_min, 1)
    if (y_min < -0.2) {
      y_breaks <- c(
        sort(seq(-0.2, y_min, by = -0.2)),
        seq(0, 1, length = 6)
      )
    } else {
      y_breaks <- c(-0.2, seq(0, 1, length = 6))
    }
  } else {
    y_lims <- c(0, 1)
    y_breaks <- seq(0, 1, length = 6)
  }
  x_breaks <- unique(data_frame$Factor)
  p <- ggplot(
    subset(data_frame, Choice == cluster),
    aes(
      x = Factor, y = Mean, group = Measure,
      color = Measure, linetype = Measure
    )
  ) +
    geom_line(stat = "identity") +
    geom_point(show.legend = FALSE) +
    geom_errorbar(aes(ymin = Mean - S.d., ymax = Mean + S.d.),
      width = .2,
      position = position_dodge(0.01), show.legend = FALSE, linetype = 1
    ) +
    expand_limits(y = y_lims) +
    scale_color_viridis_d() +
    scale_y_continuous(breaks = y_breaks) +
    scale_x_discrete(breaks = x_breaks) +
    labs(
      title = plot_title,
      subtitle = sub_title,
      x = TeX("$\\phi$"),
      y = "Measure"
    )
  suppressMessages(ggsave(filename, plot = p, compress = FALSE, device = "pdf"))
  return(p)
}

make_plot_line_issvd <- function(data_frame, cluster,
                                 x_title, y_title, filename, measures) {
  data_frame[["Measure"]] <- factor(data_frame[["Measure"]],
    labels = measures
  )

  if (min(data_frame$Mean) < 0) {
    y_min <- min(data_frame$Mean)
    y_lims <- c(y_min, 1)
    if (y_min < -0.2) {
      y_breaks <- c(
        sort(seq(-0.2, y_min, by = -0.2)),
        seq(0, 1, length = 6)
      )
    } else {
      y_breaks <- c(-0.2, seq(0, 1, length = 6))
    }
  } else {
    y_lims <- c(0, 1)
    y_breaks <- seq(0, 1, length = 6)
  }
  x_breaks <- unique(data_frame$Factor)
  p <- ggplot(
    subset(data_frame, Choice == cluster),
    aes(
      x = Factor, y = Mean, group = Measure,
      color = Measure, linetype = Measure
    )
  ) +
    geom_line(stat = "identity") +
    geom_point(show.legend = FALSE) +
    geom_errorbar(aes(ymin = Mean - S.d., ymax = Mean + S.d.),
      width = .2,
      position = position_dodge(0.01), show.legend = FALSE, linetype = 1
    ) +
    expand_limits(y = y_lims) +
    scale_color_viridis_d() +
    scale_y_continuous(breaks = y_breaks) +
    scale_x_discrete(breaks = x_breaks) +
    labs(
      x = x_title,
      y = "Measure"
    )
  p <- p +
    theme_minimal() +
    theme(
      text = element_text(size = 11), legend.position = "bottom",
      axis.text.x = element_text(angle = -45, hjust = 0.1),
      axis.ticks = element_line(size = 0.3),
      axis.line = element_line(size = 0.3)
    )
  suppressMessages(ggsave(filename, plot = p, compress = FALSE, device = "pdf"))
  return(p)
}

make_phi_agg <- function(data_frame, cluster, plot_title,
                         sub_title, x_title, y_title,
                         filename, measures, broken = TRUE, x_even = TRUE) {
  data_frame[["Measure"]] <- factor(data_frame[["Measure"]],
    labels = measures
  )

  if (broken) {
    y_lims <- c(0.6, 1)
    y_breaks <- seq(0.6, 1, length = 5)
  } else {
    y_lims <- c(0, 1)
    y_breaks <- seq(0, 1, length = 6)
  }
  if (x_even) {
    x_breaks <- unique(data_frame$Factor)
    data_frame$Factor <- factor(data_frame$Factor, labels = x_breaks)
  } else {
    x_breaks <- seq(0, 2000, 200)
  }
  p <- ggplot(
    subset(data_frame, Choice == cluster),
    aes(x = Factor, y = Mean, group = Measure, color = Measure)
  ) +
    geom_line(stat = "identity") +
    expand_limits(y = y_lims) +
    scale_color_viridis_d() +
    scale_y_continuous(breaks = y_breaks)
  if (x_even) {
    p <- p + scale_x_discrete(breaks = x_breaks)
  } else {
    p <- p + scale_x_continuous(breaks = x_breaks)
  }
  p <- p +
    labs(
      title = plot_title,
      x = TeX("$\\phi$"),
      y = "Measure"
    ) +
    theme_minimal()
  if (x_even) {
    p <- p + theme(
      text = element_text(size = 11), legend.position = "bottom",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = -45, hjust = 0.1),
      axis.ticks = element_line(size = 0.3),
      axis.line = element_line(size = 0.3)
    )
  } else {
    p <- p + theme(
      text = element_text(size = 11), legend.position = "bottom",
      legend.title = element_blank(),
      axis.ticks = element_line(size = 0.3),
      axis.line = element_line(size = 0.3)
    )
  }


  suppressMessages(ggsave(filename, plot = p, compress = FALSE, device = "pdf"))
  return(p)
}

all_plots <- function(plots, filename, n_col_plots, n_row_plots) {
  all_plots <- ggarrange(
    plots[[1]]$a, plots[[2]]$a, plots[[3]]$a, plots[[5]]$a,
    labels = c("A", "B", "C", "D"),
    font.label = list(size = 8),
    ncol = n_col_plots, nrow = n_row_plots,
    common.legend = TRUE, legend = "bottom"
  )
  suppressMessages(ggsave(filename,
    plot = all_plots,
    compress = FALSE, device = "pdf"
  ))
  return(all_plots)
}

sd_brackets <- function(mean, sd, n, measure = "other") {
  return(paste0(paste0(mean, " ("), paste0(sd, ")")))
}
collapse_rows_df <- function(df, variable) {
  group_var <- enquo(variable)
  df %>%
    group_by(!!group_var) %>%
    mutate(groupRow = seq_len(len(n))) %>%
    ungroup() %>%
    mutate(!!quo_name(group_var) := ifelse(groupRow == 1,
      as.character(!!group_var), ""
    )) %>%
    select(-c(groupRow))
}

max_index <- function(col) {
  return(which.max(as.numeric(sub(" .*", "", as.vector(col)))))
}
min_index <- function(col) {
  return(which.min(as.numeric(sub(" .*", "", as.vector(col)))))
}

all_table <- function(df, measures, filename, repeats, col_names_list,
                      subheading, table_head, ndig = 4,
                      feat = "Method", group = FALSE) {
  df <- mutate_at(df, c(feat, "Measure"), as.factor)
  df <- pivot_wider(df, names_from = Measure, values_from = c("Mean", "S.d."))
  # round to 4 digits
  df <- mutate_if(
    df, is.numeric,
    function(x) format(round(x, ndig), nsmall = ndig)
  )
  # add info in brackets
  mean_col_names <- paste0("Mean_", measures)
  sd_col_names <- paste0("S.d._", measures)
  # name of mean column
  for (i in seq_len(length(measures))) {
    df[[measures[i]]] <- sd_brackets(
      df[[mean_col_names[i]]],
      df[[sd_col_names[i]]], repeats, measures[i]
    )
  }

  df_all <- subset(df, select = c("Method", "Choice", measures))
  df_overall <- subset(df, select = c(feat, measures))
  colnames(df_all) <- c(feat, "Choice", measures)
  colnames(df_overall) <- col_names_list
  df_overall <- as.data.frame(t(df_overall))
  df_overall <- df_overall[-1, ]
  colnames(df_overall) <- unique(df[[feat]])
  n <- ncol(df_overall)
  vec <- c(" " = 1)
  vec[[table_head]] <- n
  df_overall2 <- kable(df_overall)
  df_overall2 <- kbl(df_overall, booktabs = T, "latex", escape = FALSE)
  df_overall2 <- add_header_above(df_overall2, header = vec)
  df_overall2 <- kable_styling(add_footnote(df_overall2, subheading, escape = FALSE))
  # save
  sink(paste0(filename, "_overall.txt"))
  print(df_overall2)
  sink()

  # collapse method rows
  if (group) {
    group_tab <- kable_styling(kbl(df_overall, booktabs = T, "latex"))
    # select groups
    feats <- unique(df_overall[[feat]])
    for (i in seq_len(length(feats))) {
      group_tab <- pack_rows(group_tab, feats[i], 3 * i - 2, 3 * i)
    }
    sink(paste0(filename, "_all.txt"))
    print(group_tab)
    sink()
  }
  return(df_overall2)
}



measure_table <- function(df, measure_val, filename, repeats, col_names, subheading, table_head, ndig = 4) {
  factors <- unique(df$Factor)
  df <- subset(df, Measure == measure_val)
  df <- mutate_at(df, c("Measure"), as.factor)
  df <- pivot_wider(df, names_from = Factor, values_from = c("Mean", "S.d."))
  # round to 4 digits
  df <- mutate_if(df, is.numeric, function(x) format(round(x, ndig), nsmall = ndig))
  # add info in brackets
  mean_col_names <- paste0("Mean_", factors)
  sd_col_names <- paste0("S.d._", factors)
  # name of mean column
  for (i in seq_len(length(factors))) {
    df[[as.character(factors[i])]] <- sd_brackets(
      df[[mean_col_names[i]]],
      df[[sd_col_names[i]]], repeats, factors[i]
    )
  }
  c_factors <- as.character(factors)
  df_all <- subset(df, select = c("Method", c_factors))
  df_overall <- subset(df, select = c("Method", c_factors))
  colnames(df_all) <- c("Method", col_names)
  colnames(df_overall) <- c("", col_names)
  # subset only overall results and remove extras
  # make the largest values down a column bold
  for (i in seq_len(length(factors))) {
    max_id <- max_index(df_overall[[c_factors[i]]])
    cond <- seq_len(length(df_overall[, 1][[1]])) == max_id
    df_overall[[c_factors[i]]] <- text_spec(df_overall[[c_factors[i]]], "latex",
      bold = cond
    )
  }
  n <- ncol(df_overall) - 1
  vec <- c(" " = 1)
  vec[[table_head]] <- n
  df_overall2 <- kable(df_overall)
  df_overall2 <- kbl(df_overall, booktabs = T, "latex", escape = FALSE)
  df_overall2 <- add_header_above(df_overall2, header = vec)
  df_overall2 <- kable_styling(add_footnote(df_overall2, subheading, escape = FALSE))
  # save
  sink(paste0(filename, "_overall.txt"))
  print(df_overall2)
  sink()
  return(list("table" = df_all, "latex" = df_overall2))
}

combine_tables <- function(tables, filename, measures, table_head, subheading, col_names_tables) {
  full_table <- rbind.data.frame(
    tables[[1]]$table, tables[[2]]$table,
    tables[[3]]$table, tables[[4]]$table, tables[[5]]$table
  )
  # collapse method rows
  n <- ncol(full_table) - 1
  n2 <- length(unique(full_table[["Method"]]))
  colnames(full_table) <- c("", col_names_tables)
  group_tab <- kable_styling(kbl(full_table, booktabs = T, "latex"))
  vec <- c(" " = 1)
  vec[[table_head]] <- n
  group_tab <- add_header_above(group_tab, header = vec)
  group_tab <- kable_styling(add_footnote(group_tab, subheading, escape = FALSE))
  # select groups
  for (i in seq_len(length(measures))) {
    group_tab <- pack_rows(group_tab, measures[i], n2 * i - 3, n2 * i)
  }

  sink(paste0(filename, "_all.txt"))
  print(group_tab)
  sink()
  return(group_tab)
}

corr_table <- function(results, filename) {
  measures <- c("F_score", "Rel", "Rec", "BiS", "CSR")
  corr_tab <- matrix(0, nrow = 3, ncol = 2)
  for (i in 1:3) {
    for (j in 1:2) {
      vec1 <- results$Mean[results$Measure == measures[i]]
      vec2 <- results$Mean[results$Measure == measures[3 + j]]
      corr_tab[i, j] <- cor(vec1, vec2)
    }
  }
  rownames(corr_tab) <- c("F score", "Relevance", "Recovery")
  colnames(corr_tab) <- c("BiS", "CSR")
  latex_tab <- kbl(corr_tab, booktabs = T, "latex", escape = FALSE)
  # latex_tab <- kable_styling(add_footnote(latex_tab, subheading))
  # save
  sink(paste0(filename, "corr_tab.txt"))
  print(latex_tab)
  sink()
  return(latex_tab)
}
