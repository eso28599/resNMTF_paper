## functions to remove spurious biclusters
get_thresholds <- function(data, output_f, repeats) {
  # Finds the threshold for removal of spurious biclusters.
  # Args:
  #  data: list of matrices
  #  output_f: list of matrices
  #  repeats: minimum value of 2
  # Returns:
  #  list of thresholds
  n_views <- length(data)
  n_clusts <- dim(output_f[[1]])[2]
  k_input <- n_clusts * rep(1, length = n_views)
  x_mess <- vector(mode = "list", length = repeats)
  for (n in 1:repeats) {
    data_messed <- vector(mode = "list", length = n_views)
    for (i in 1:n_views) {
      # correct shuffling
      dims <- dim(data[[i]])
      data_messed[[i]] <- matrix(
        sample(data[[i]]),
        dims[1], dims[2]
      )
      while (any(colSums(data_messed[[i]]) == 0) ||
        any(rowSums(data_messed[[i]]) == 0)) {
        data_messed[[i]] <- matrix(
          sample(data[[i]]),
          dims[1], dims[2]
        )
      }
    }
    results <- restMultiNMTF_run(
      data = data_messed,
      k_vec = k_input, no_clusts = TRUE, stability = FALSE
    )
    x_mess[[n]] <- results$output_f
  }
  avg_score <- c()
  max_score <- c()
  data <- vector(mode = "list", length = n_views)
  dens_list <- vector(mode = "list", length = n_views)
  d <- 1
  for (i in 1:n_views) {
    scores <- c()
    for (j in 1:max(repeats - 1, 1)) {
      data[[i]] <- cbind(data[[i]], x_mess[[j]][[i]])
      for (k in 1:n_clusts) {
        # jth repeat, ith view, kth cluster
        x1 <- ((x_mess[[j]])[[i]])[, k]
        for (l in (j + 1):repeats) {
          for (m in 1:n_clusts) {
            x2 <- ((x_mess[[l]])[[i]])[, m]
            max_val <- max(x1, x2)
            d1 <- density(x1, from = 0, to = max_val)
            d2 <- density(x2, from = 0, to = max_val)
            d1$y[d1$x > max(x1)] <- 0
            d2$y[d2$x > max(x2)] <- 0
            dens_list[[d]] <- d1
            scores <- c(
              scores,
              suppressMessages(JSD(rbind(d1$y, d2$y),
                unit = "log2", est.prob = "empirical"
              ))
            )
          }
        }
      }
    }
    data[[i]] <- cbind(data[[i]], x_mess[[repeats]][[i]])
    avg_score <- c(avg_score, mean(scores))
    dens <- density(scores)
    max_score <- c(max_score, dens$x[which.max(dens$y)])
  }
  return(list(
    "avg_score" = avg_score,
    "max_score" = max_score,
    "scores" = scores,
    "data" = data,
    "dens" = dens_list
  ))
}

check_biclusters <- function(data, output_f, repeats) {
  # calculate JSD between returned F and noise
  # returns scores, avg score and max score
  n_views <- length(data)
  n_clusts <- dim(output_f[[1]])[2]
  # updated results
  scores <- matrix(0, nrow = n_views, ncol = n_clusts)
  thresholds <- get_thresholds(data, output_f, repeats)
  for (i in 1:n_views) {
    x_noise <- thresholds$data[[i]]
    for (k in 1:n_clusts) {
      x <- output_f[[i]][, k]
      scores[i, k] <- mean(apply(
        x_noise, 2,
        function(y) jsd_calc(x, y)
      ))
    }
  }
  return(list(
    "score" = scores,
    "avg_threshold" = thresholds$avg_score,
    "max_threshold" = thresholds$max_score
  ))
}

obtain_biclusters <- function(data, output_f,
                              output_g, output_s, repeats, distance) {
  # Obtains biclustering from ResNMTF factorisation
  # Args:
  #  data: list of data matrices
  #  output_f: list of F matrices
  #  output_g: list of G matrices
  #  output_s: list of S matrices
  #  repeats: minimum value of 2
  #  distance: distance metric to use, can be . Default is euclidean.
  # Returns:
  #  list containing row and column clusterings, and the bisilhouette score of the biclustering

  n_views <- length(output_f)
  row_clustering <- vector("list", length = n_views)
  col_clustering <- vector("list", length = n_views)

  # assign biclusters
  biclusts <- check_biclusters(data, output_f, repeats)
  for (i in 1:n_views) {
    row_clustering[[i]] <- apply(
      output_f[[i]],
      2, function(x) as.numeric(x > (1 / dim(output_f[[i]])[1]))
    )

    col_clustering[[i]] <- apply(
      output_g[[i]],
      2, function(x) as.numeric(x > (1 / dim(output_g[[i]])[1]))
    )
  }
  bisil <- c()
  # update realtions and
  # set biclusters that aren't strong enough to 0
  # and if bicluster is empty set row and cols to 0
  for (i in 1:n_views) {
    indices <- (((biclusts$score[i, ]) < biclusts$max_threshold[i]) |
      ((biclusts$score[i, ]) == 0))
    relations <- apply(output_s[[i]], 1, which.max)
    new_indices <- indices[relations] # i==0, col cluster i isn't a bicluster
    row_clustering[[i]] <- row_clustering[[i]][, relations]
    row_clustering[[i]][, new_indices] <- 0
    col_clustering[[i]][, new_indices] <- 0
    bisil <- c(
      bisil,
      bisilhouette(data[[i]], row_clustering[[i]], col_clustering[[i]], method = distance)$bisil
    )
  }
  # calculate overall bisil
  bisil <- ifelse(sum(bisil) == 0, 0, mean(bisil[bisil != 0]))
  return(list(
    "row_clustering" = row_clustering,
    "col_clustering" = col_clustering, "bisil" = bisil
  ))
}
