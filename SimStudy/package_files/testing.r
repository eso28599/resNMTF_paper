source("/Users/ellaorme/Documents/resNMTF/SimStudy/package_files/main.r")
source("/Users/ellaorme/Documents/resNMTF/SimStudy/package_files/obtain_bicl.r")
source("/Users/ellaorme/Documents/resNMTF/SimStudy/package_files/stability_analysis.r")
source("/Users/ellaorme/Documents/resNMTF/SimStudy/package_files/update_steps.r")
source("/Users/ellaorme/Documents/resNMTF/SimStudy/package_files/utils.r")

row_clusters <- cbind(
  rbinom(100, 1, 0.5),
  rbinom(100, 1, 0.5),
  rbinom(100, 1, 0.5)
)
col_clusters <- cbind(
  rbinom(50, 1, 0.4),
  rbinom(50, 1, 0.4),
  rbinom(50, 1, 0.4)
)
n_col <- 50
data_1 <- list(
  row_clusters %*% diag(c(5, 5, 5)) %*% t(col_clusters) +
    abs(matrix(rnorm(100 * n_col), 100, n_col)),
  row_clusters %*% diag(c(5, 5, 5)) %*% t(col_clusters) +
    abs(0.01 * matrix(rnorm(100 * n_col), 100, n_col))
)


start.time <- Sys.time()
res <- restMultiNMTF_run(data_1, k_max = 4)
end.time <- Sys.time()
end.time - start.time

library(bisilhouette)

relevance_results(
  res$row_clusters[[1]],
  res$col_clusters[[1]],
  row_clusters, col_clusters
)
relevance_results(
  res$row_clusters[[2]],
  res$col_clusters[[2]],
  row_clusters, col_clusters
)
jaccard_res(
  res$row_clusters[[1]],
  res$col_clusters[[1]],
  row_clusters, col_clusters
)
jaccard_res(
  res$row_clusters[[2]],
  res$col_clusters[[2]],
  row_clusters, col_clusters
)

source("SimStudy/RunSim/Functions/data_generation.r")

data <- multi_view(rep(200, 3), rep(200, 3), 3, 5, 5, row_e = 1, col_e = 1, row_o = 0, col_o = 0, row_same_shuffle = TRUE, col_same_shuffle = FALSE, seed = FALSE, file_path = NA)
res_data <- restMultiNMTF_run(data$data_views, k_max = 4)
jaccard_res(
  res_data$row_clusters[[1]],
  res_data$col_clusters[[1]],
  data$truth_rows[[1]], data$truth_cols[[1]]
)

jaccard_main <- function(samps, feats, row_c, col_c, true_r, true_c, m, n) {
  # initialise storage of jaccard index between pairs
  jac_mat <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    r_i <- samps[row_c[, i] == 1]
    c_i <- feats[col_c[, i] == 1]
    m_i <- cart_prod(r_i, c_i)
    for (j in 1:n) {
      tr_i <- samps[true_r[, j] == 1]
      tc_i <- feats[true_c[, j] == 1]
      m_j <- cart_prod(tr_i, tc_i)
      jac_mat[i, j] <- jaccard_func(m_i, m_j)
    }
  }
  return(jac_mat)
}

relevance_results <- function(row_c, col_c, true_r, true_c) {
  m <- ncol(row_c)
  n <- ncol(true_r)
  m_0 <- sum(colSums(row_c) != 0) # no of clusters actually detected
  n_0 <- sum(colSums(true_r) != 0) # no of true clusters

  # edge cases
  # if either:
  #     - no biclusters detected but some are present
  #     - no biclusters present but some are detected
  # return 0
  if ((m_0 == 0 && n_0 != 0) || (n_0 == 0 && m_0 != 0)) {
    return(0)
  }
  # if no biclusters present and none detected
  # return 1
  if (m_0 == 0 && n_0 == 0) {
    return(1)
  }
  samps <- seq_along(row_c[, 1])
  feats <- seq_along(col_c[, 1])
  # initialise storage of jaccard index between pairs
  jac_mat <- jaccard_main(samps, feats, row_c, col_c, true_r, true_c, m, n)
  return(apply(jac_mat, 2, max))
}
jaccard_func <- function(a, b) {
  # calculate jaccard between two vectors
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  if (union == 0) {
    return(0)
  } else {
    return(intersection / union)
  }
}
