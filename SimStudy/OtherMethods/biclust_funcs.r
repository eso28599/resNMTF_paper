# library(biclust)
library(bisilhouette)
apply_biclust_method <- function(data, method, name,
                                 transposed = FALSE, number_clusts = 6) {
  row_clusters <- vector("list", length = length(data))
  col_clusters <- vector("list", length = length(data))
  for (i in seq_along(data)) {
    res <- suppressMessages(biclust::biclust(data[[i]], method = method))
    if (anyNA(res@RowxNumber)) {
      row_clusters[[i]] <- matrix(0, nrow = nrow(data[[i]]), ncol = 1)
      col_clusters[[i]] <- matrix(0, nrow = ncol(data[[i]]), ncol = 1)
    } else {
      rows <- apply(res@RowxNumber, 2, as.numeric)
      cols <- apply(res@NumberxCol, 2, as.numeric)
      if (!is.matrix(cols)) {
        col_clusters[[i]] <- matrix(cols, nrow = length(cols), ncol = 1)
        row_clusters[[i]] <- matrix(rows, nrow = length(rows), ncol = 1)
      } else {
        col_clusters[[i]] <- t(cols)
        row_clusters[[i]] <- rows
      }
    }
  }
  if (transposed) {
    return(list(row_clusters = col_clusters, col_clusters = row_clusters))
  }
  return(list(row_clusters = row_clusters, col_clusters = col_clusters))
}
biclust_results <- function(data, labels, method, name,
                            transposed = FALSE, repeats = 5,
                            number_clusts = 6) {
  n_cols <- (length(data) + 1) * 6
  res_mat <- matrix(0, nrow = repeats, ncol = n_cols)
  for (i in seq_len(repeats)) {
    results <- apply_biclust_method(
      data, method, name,
      transposed, number_clusts
    )
    jaccs <- all_jaccs(labels, results$row_clusters)
    n_views <- length(data)
    res_mat[i, ] <- c(jaccs, rep(0, (n_views + 1) * (3)))
  }
  res_mat <- cbind(1:5, rep(name, repeats), res_mat)
  return(res_mat)
}
