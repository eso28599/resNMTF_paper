source("SimStudy/Functions/data_generation.r")
source("SimStudy/Functions/extra_funcs.r")
library(bisilhouette)
library(resnmtf)
# generate data
set.seed(123)
n_views <- 3
row_cl_dims <- rep(200, n_views)
# 100, 50,250 features respectively
col_cl_dims <- c(100, 50, 250)
save_data(row_cl_dims, col_cl_dims, 5,
  "Exploration/visual_data", 5,
  col_same_shuffle = FALSE
)

# import data
data <- import_matrix_og("Exploration/visual_data/data.xlsx")
true_rows <- import_matrix_og("Exploration/visual_data/true_rows.xlsx")
true_cols <- import_matrix_og("Exploration/visual_data/true_cols.xlsx")

# set seed
rows_shuffled <- cbind(
  true_rows[[1]][, c(1, 2, 3)],
  sample(true_rows[[1]][, 4]),
  sample(true_rows[[1]][, 5])
)

# true plot
bisil_plot(data[[1]], true_rows[[1]], true_cols[[1]],
  filename = "Exploration/visual_data/true_bisil_plot.pdf",
  method = "euclidean",
  w = 7, h = 7
)
# plot for shuffled cols for two biclusters
bisil_plot(data[[1]], rows_shuffled, true_cols[[1]],
  filename = "Exploration/visual_data/shuffled_bisil_plot.pdf",
  method = "euclidean",
  w = 7, h = 7
)

phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3)] <- 200
phi[2, c(3)] <- 200
results <- apply_resnmtf(data,
  phi = phi, k_min = 5, k_max = 6,
  stability = FALSE
)
write.csv(results$All_Error,
  "Exploration/visual_data/resnmtf_error.csv",
  row.names = FALSE
)
