source("SimStudy/Functions/data_generation.r")
source("SimStudy/Functions/extra_funcs.r")
library(bisilhouette)
library(resnmtf)
library(ggplot2)

# generate data
set.seed(100)
n_views <- 3
row_cl_dims <- rep(200, n_views)
# 100, 50,250 features respectively
col_cl_dims <- c(100, 50, 250)
phi_mat <- matrix(0, nrow = n_views, ncol = n_views)
phi_mat[1, c(2, 3)] <- 200
phi_mat[2, c(3)] <- 200
for (k in 3:6) {
  save_data(row_cl_dims, col_cl_dims, k,
    paste0("Exploration/convergence_data/", k),
    5,
    col_same_shuffle = FALSE
  )
  data <- import_matrix_og(
    paste0("Exploration/convergence_data/", k, "/data.xlsx")
  )
  results <- apply_resnmtf(data,
    phi = phi_mat,
    k_vec = rep(k, n_views),
    stability = FALSE
  )$All_Error
  write.csv(results,
    paste0("Exploration/convergence_data/", k, "/all_results.csv"),
    row.names = FALSE
  )
}

errors_3 <- read.csv("Exploration/convergence_data/3/all_results.csv")
errors_4 <- read.csv("Exploration/convergence_data/4/all_results.csv")
errors_5 <- read.csv("Exploration/convergence_data/5/all_results.csv")
errors_6 <- read.csv("Exploration/convergence_data/6/all_results.csv")

n_its <- max(
  nrow(errors_3),
  nrow(errors_4),
  nrow(errors_5),
  nrow(errors_6)
)
all_errors <- matrix(0, nrow = n_its, ncol = 4)
all_errors[seq_len(nrow(errors_3)), 1] <- errors_3[, 1]
all_errors[seq_len(nrow(errors_4)), 2] <- errors_4[, 1]
all_errors[seq_len(nrow(errors_5)), 3] <- errors_5[, 1]
all_errors[seq_len(nrow(errors_6)), 4] <- errors_6[, 1]
colnames(all_errors) <- c("3", "4", "5", "6")
write.csv(all_errors,
  "Exploration/convergence_data/all_errors.csv",
  row.names = FALSE
)


all_errors <- read.csv("Exploration/convergence_data/all_errors.csv")
n_its <- nrow(all_errors)
df <- data.frame("Iteration" = 1:n_its, all_errors)
colnames(df) <- c("Iteration", "3", "4", "5", "6")
df_long <- reshape2::melt(df,
  id.vars = "Iteration",
  variable.name = "K", value.name = "Error"
)
df_long <- subset(df_long, all_errors != 0)
p <- ggplot(data = df_long, aes(x = Iteration, y = Error, group = K)) +
  geom_line(aes(color = K), linewidth = 0.5) +
  labs(x = "Iteration", y = "Error") +
  scale_color_viridis_d(name = "K") +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_text(size = 10))

ggsave("Exploration/convergence_data/error_plot.pdf",
  plot = p, compress = FALSE, device = "pdf", width = 5, height = 4
)

# import plot for 3cancers
cancers_errors <- read.csv("RealData/3cancers/data/stability/errors_1.csv")
df_cancers <- data.frame(
  "Iteration" = seq_len(nrow(cancers_errors)),
  "Error" = cancers_errors[, 1]
)
p_cancers <- ggplot(data = df_cancers, aes(x = Iteration, y = Error)) +
  geom_line(color = "blue", linewidth = 0.5) +
  labs(x = "Iteration", y = "Error") +
  theme_minimal()

ggsave("RealData/3cancers/data/stability/error_plot.pdf",
  plot = p_cancers, compress = FALSE, device = "pdf", width = 5, height = 4
)
