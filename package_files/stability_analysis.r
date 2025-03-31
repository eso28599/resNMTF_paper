# ---------------------------------
# functions for stability selection
# ---------------------------------

jaccard_results <- function(row_c, col_c, true_r, true_c) {
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
  samps <- seq_along(row_c)
  feats <- seq_along(col_c)
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
      jac_mat[i, j] <- jaccard(m_i, m_j)
    }
  }
  return(apply(jac_mat, 2, max))
}

#' Test whether there are any columns/rows with only zeros
test_cond <- function(data, attempt) {
  if (attempt == 1) {
    return(TRUE)
  }
  return(any(unlist(lapply(
    data,
    function(x) {
      any(colSums(x) == 0) | any(rowSums(x) == 0)
    }
  ))))
}


stability_check <- function(data, init_s, results,
                            k, phi, xi, psi, n_iters,
                            repeats, no_clusts, distance, sample_rate = 0.9,
                            n_stability = 5, stab_thres = 0.6,
                            stab_test = FALSE) {
  # check whether stability check even needs to be done
  # no_clusts_detected
  n_c <- sum(as.numeric(lapply(
    results$row_clusters,
    function(x) sum(colSums(x))
  )))
  if (n_c == 0) {
    print("No biclusters detected!")
    return(results)
  }
  n_views <- length(data)
  dim <- dim(data[[1]])
  # initialise storage of results
  jacc <- matrix(0, nrow = n_views, ncol = k)
  for (t in 1:n_stability) {
    new_data <- vector(mode = "list", length = n_views)
    row_samples <- vector(mode = "list", length = n_views)
    col_samples <- vector(mode = "list", length = n_views)
    # turn this into a function to be used with lapply
    attempt <- 1
    # need to check that
    while (test_cond(new_data, attempt)) {
      if (attempt == 20) {
        print("Unable to perform stability analysis due to sparsity of data.")
        return(results)
      }
      row_samples[[1]] <- sample(dim[1], (dim[1] * sample_rate))
      col_samples[[1]] <- sample(dim[2], (dim[2] * sample_rate))
      new_data[[1]] <- data[[1]][row_samples[[1]], col_samples[[1]]]
      if (any(colSums(new_data[[1]]) == 0) | any(rowSums(new_data[[1]]) == 0)) {
        zeros_cols <- colSums(new_data[[1]]) != 0
        zeros_rows <- rowSums(new_data[[1]]) != 0
        row_samples[[1]] <- row_samples[[1]][zeros_rows]
        col_samples[[1]] <- col_samples[[1]][zeros_cols]
        new_data[[1]] <- data[[1]][row_samples[[1]], col_samples[[1]]]
      }
      if (n_views > 1) {
        for (i in 2:n_views) {
          dims <- dim(data[[i]])
          if ((dims[1]) == dim[1]) {
            row_samples[[i]] <- row_samples[[1]]
          } else {
            row_samples[[i]] <- sample(dims[1], (dims[1] * sample_rate))
          }
          if ((dims[2]) == dim[2]) {
            col_samples[[i]] <- col_samples[[1]]
          } else {
            col_samples[[i]] <- sample(dims[2], (dims[2] * sample_rate))
          }
          new_data[[i]] <- data[[i]][row_samples[[i]], col_samples[[i]]]
          if (any(colSums(new_data[[i]]) == 0) ||
            any(rowSums(new_data[[i]]) == 0)) {
            zeros_cols <- colSums(new_data[[i]]) != 0
            zeros_rows <- rowSums(new_data[[i]]) != 0
            if ((dims[1]) == dim[1]) {
              for (p in 1:i) {
                row_samples[[p]] <- row_samples[[p]][zeros_rows]
              }
            } else {
              row_samples[[i]] <- row_samples[[i]][zeros_rows]
            }
            if ((dims[2]) == dim[2]) {
              for (p in 1:i) {
                col_samples[[p]] <- col_samples[[p]][zeros_cols]
              }
            } else {
              col_samples[[i]] <- col_samples[[i]][zeros_cols]
            }
            for (p in 1:i) {
              new_data[[p]] <- data[[p]][row_samples[[p]], col_samples[[p]]]
            }
          }
        }
      }
      attempt <- attempt + 1
    }
    new_results <- rest_multi_nmtf_inner(
      new_data,
      k_vec = k,
      phi = phi,
      xi = xi,
      psi = psi,
      n_iters = n_iters,
      repeats = repeats,
      distance = distance
    )
    # extract results
    for (i in 1:n_views) {
      jacc[i, ] <- jacc[i, ] + jaccard_results(
        new_results$row_clusters[[i]],
        new_results$col_clusters[[i]],
        results$row_clusters[[i]][row_samples[[i]], ],
        results$col_clusters[[i]][col_samples[[i]], ]
      )
    }
  }
  jacc <- jacc / n_stability
  if (stab_test) {
    return(list("res" = results, "jacc" = jacc))
  } else {
    for (i in 1:n_views) {
      # set clusters not deemed stable to have 0 members
      results$row_clusters[[i]][, jacc[i, ] < stab_thres] <- 0
      results$col_clusters[[i]][, jacc[i, ] < stab_thres] <- 0
    }
    return(results)
  }
}
