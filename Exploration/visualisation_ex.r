source("SimStudy/Functions/data_generation.r")
source("SimStudy/Functions/extra_funcs.r")
library(bisilhouette)
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
data <- import_matrix("Exploration/visual_data/data.xlsx")
true_rows <- import_matrix("Exploration/visual_data/true_rows.xlsx")
true_cols <- import_matrix("Exploration/visual_data/true_cols.xlsx")

# set seed
rows_shuffled <- cbind(
  true_rows[[1]][, c(1, 2, 3)],
  sample(true_rows[[1]][, 4]),
  sample(true_rows[[1]][, 5])
)

# true plot
bisil_plot(data[[1]], true_rows[[1]], true_cols[[1]],
  filename = "Exploration/visual_data/true_bisil_plot.pdf"
)
# plot for shuffled cols for two biclusters
bisil_plot(data[[1]], rows_shuffled, true_cols[[1]],
  filename = "Exploration/visual_data/shuffled_bisil_plot.pdf"
)


# overlap example plot
# path <- "investigation/overlap_nonex_example"
# data <- multi_view(row_cl_dims, col_cl_dims, 5, 5, 5, row_e = 0.9, col_e =0.9, row_o = 0.2, col_o = 0.2, row_same_shuffle=TRUE, col_same_shuffle=FALSE, seed=123, file_path=path)
