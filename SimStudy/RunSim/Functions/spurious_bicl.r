# removal of spurious biclusters
# source("utils.r")
# Libraries
library(philentropy)
library("Matrix")
library("clue")
library("aricode")
library("rio")
library("eList")
library(foreach)
library(doParallel)
library(doSNOW)
library(MASS)
## functions to remove spurious biclusters
### functions for calculating JSD
jsd_calc <- function(x1, x2) {
  # calculate the Jensen Shannon divergence between
  # vectors x1 and x2
  max_val <- max(x1, x2)
  d1 <- density(x1, from = 0, to = max_val)
  d2 <- density(x2, from = 0, to = max_val)
  d1$y[d1$x > max(x1)] <- 0
  d2$y[d2$x > max(x2)] <- 0
  return(suppressMessages(JSD(rbind(d1$y, d2$y),
    unit = "log2",
    est.prob = "empirical"
  )))
}


get_thresholds <- function(Xinput, Foutput, repeats) {
  # Finds the threshold for removal of spurious biclusters.
  # Args:
  #  Xinput: list of matrices
  #  Foutput: list of matrices
  #  repeats: minimum value of 2
  # Returns:
  #  list of thresholds
  n_views <- length(Xinput)
  n_clusts <- dim(Foutput[[1]])[2]
  k_input <- n_clusts * rep(1, length = n_views)
  x_mess <- vector(mode = "list", length = repeats)
  for (n in 1:repeats) {
    data_messed <- vector(mode = "list", length = n_views)
    for (i in 1:n_views) {
      # correct shuffling
      dims <- dim(Xinput[[i]])
      data_messed[[i]] <- matrix(
        sample(Xinput[[i]]),
        dims[1], dims[2]
      )
      while (any(colSums(data_messed[[i]]) == 0) |
        any(rowSums(data_messed[[i]]) == 0)) {
        data_messed[[i]] <- matrix(
          sample(Xinput[[i]]),
          dims[1], dims[2]
        )
      }
    }
    results <- restMultiNMTF_run(
      Xinput = data_messed,
      KK = k_input, no_clusts = TRUE, stability = FALSE
    )
    x_mess[[n]] <- results$Foutput
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
                                   unit = "log2", est.prob = "empirical"))
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

check_biclusters <- function(Xinput, Foutput, repeats) {
  # calculate JSD between returned F and noise
  # returns scores, avg score and max score
  n_views <- length(Xinput)
  n_clusts <- dim(Foutput[[1]])[2]
  # updated results
  scores <- matrix(0, nrow = n_views, ncol = n_clusts)
  thresholds <- get_thresholds(Xinput, Foutput, repeats)
  for (i in 1:n_views) {
    x_noise <- thresholds$data[[i]]
    for (k in 1:n_clusts) {
      x <- Foutput[[i]][, k]
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
