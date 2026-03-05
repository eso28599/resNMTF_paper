library(openxlsx)
library(MASS)
library(Matrix)
source("SimStudy/Functions/data_generation.r")

# original overlap used
data <- multi_view(c(200), c(200), 5, 5, 5,
  row_e = 0.8, col_e = 0.8, row_o = 0.3, col_o = 0.3,
  row_same_shuffle = TRUE, col_same_shuffle = TRUE,
  seed = 1, file_path = "Exploration/visual_data/original_overlap"
)

# additive overlap
one_view_adv <- function(
    row_dims, col_dims, k, noise, signal,
    row_e = 1, col_e = 1, row_o = 0, col_o = 0) {
  #' Description:
  #'  Generate a single view of the data
  #' Arguments:
  #' row_dims: number of rows in this view
  #'  col_dims: number of columns in this view
  #'  k: number of clusters
  #'  row_e: portion of rows in biclusters
  #'  col_e: portion of columns in biclusters
  #'  row_o: portion of overlap in rows
  #'  col_o: portion of overlap in columns
  #'  signal: mean of the signal in this view
  #'  noise: variance of the noise added to views
  #' Returns:
  #'  view: matrix of this dataview
  #'  truth_row: vector indicating membership of row clusters
  #'  truth_col: vector indicating membership of col clusters

  # generate noise
  x_noise <- mvrnorm(
    n = row_dims, mu = rep(0, col_dims),
    Sigma = diag(noise, col_dims)
  )

  # Create a list as input for block-diagonal data generation
  # list length of no of clusters in this first view

  # generate a mvn with n=number of individuals of view and no of features
  # equal to number of features of the view
  # mean of 5 for each column
  # covariance matrix is identity - each feature is independent
  # define vectors indicating true row/column cluster membership for each view
  true_row <- matrix(0, nrow = row_dims, ncol = k)
  true_col <- matrix(0, nrow = col_dims, ncol = k)
  x_view <- matrix(0, nrow = row_dims, ncol = col_dims)
  dims <- get_dims(row_dims, col_dims, k, row_e, col_e, row_o, col_o)
  row_start <- dims$Row_s
  row_end <- dims$Row_e
  col_start <- dims$Col_s
  col_end <- dims$Col_e
  # row_start
  for (i in 1:k) {
    n_r <- (row_end[i] - row_start[i] + 1)
    n_c <- (col_end[i] - col_start[i] + 1)
    if (n_r != 0) {
      # original overlap used in the manuscript
      # x_view[(row_start[i]):(row_end[i]), (col_start[i]):(col_end[i])] <-
      #   mvrnorm(n = n_r, mu = rep(signal, n_c), Sigma = diag(n_c))
      # additive overlap
      x_view[(row_start[i]):(row_end[i]), (col_start[i]):(col_end[i])] <-
        x_view[(row_start[i]):(row_end[i]), (col_start[i]):(col_end[i])] +
        mvrnorm(n = n_r, mu = rep(signal, n_c), Sigma = diag(n_c))
      true_row[(row_start[i]):(row_end[i]), i] <- 1
      true_col[(col_start[i]):(col_end[i]), i] <- 1
    }
  }
  # add noise to first view
  x_view <- abs(x_view) + abs(x_noise)
  return(list(view = x_view, truth_row = true_row, truth_col = true_col))
}

# additive overlap examples
# supplementary figure
data <- multi_view(c(200), c(200), 5, 5, 5,
  row_e = 0.8, col_e = 0.8, row_o = 0.3, col_o = 0.3,
  row_same_shuffle = TRUE, col_same_shuffle = TRUE,
  seed = 1, file_path = "Exploration/visual_data/additive_overlap"
)

## plots for presentations
# generate data with different combinations of overlap and exhaustiveness

# non-overlap but exhaustive
data <- multi_view(c(200), c(600), 3, 2, 3,
  row_e = 1, col_e = 1, row_o = 0, col_o = 0,
  row_same_shuffle = TRUE, col_same_shuffle = TRUE,
  seed = 1, file_path = "Exploration/visual_data/nonoverlap_ex"
)
# non-overlap but exhaustive
data1 <- multi_view(c(200), c(600), 3, 2, 3,
  row_e = 1, col_e = 1, row_o = 0.3, col_o = 0.3,
  row_same_shuffle = TRUE, col_same_shuffle = TRUE,
  seed = 1, file_path = "Exploration/visual_data/overlap_ex"
)
# non-overlap but non-exhaustive
data1 <- multi_view(c(200), c(600), 3, 2, 3,
  row_e = 0.8, col_e = 0.8, row_o = 0, col_o = 0,
  row_same_shuffle = TRUE, col_same_shuffle = TRUE,
  seed = 1, file_path = "Exploration/visual_data/nonoverlap_nonex"
)
# overlap and non-exhaustive
data <- multi_view(c(200), c(600), 3, 2, 3,
  row_e = 0.8, col_e = 0.8, row_o = 0.3, col_o = 0.3,
  row_same_shuffle = TRUE, col_same_shuffle = TRUE,
  seed = 1, file_path = "Exploration/visual_data/overlap_nonex"
)
