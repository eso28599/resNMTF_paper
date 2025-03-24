library(philentropy)
library(Matrix)
library(clue)
library(aricode)
library(rio)
library(eList)
library(foreach)
library(doParallel)
library(doSNOW)
library(MASS)

single_alt_l1_normalisation <- function(Xmatrix) {
  # normalises a matrix so that the l1 norm of each column is 1
  Q <- diag(colSums(Xmatrix))
  # solve Q
  newMatrix <- Xmatrix %*% solve(Q)
  return(list("Q" = Q, "newMatrix" = newMatrix))
}

check <- function(matrix) {
  if (sum(colSums(matrix) != 0) > 1) {
    matrix <- matrix[, colSums(matrix) != 0]
  }
  n_clusts <- ncol(matrix)
  equal <- diag(n_clusts)

  for (i in 1:(n_clusts - 1)) {
    for (j in (i + 1):n_clusts) {
      check <- all(matrix[, i] == matrix[, j])
      equal[i, j] <- check
      equal[j, i] <- check
    }
  }
  return(nrow(unique(equal)) == 2)
}

BisilScoreInner <- function(Xinput, row_clustering, col_clustering, method="euclidean"){
  #' Calculate the bisilhouette score without repeats
  #'
  #' Args:
  #'  Xinput: data matrix, shape (N, p).
  #'  row_clustering: binary matrix indicating row clustering, shape(N, k).
  #'  col_clustering: binary matrix indicating column clustering, shape(p, k).
  #'  method: distance metric to use, str. Default is "euclidean".
  #'
  #' Returns:
  #'  bisil: bisilhouette score, float.
  
  n_clusts <- ncol(row_clustering)
  n_clusts_row <- n_clusts
  bisil_score <- rep(0, length = n_clusts)

  # ensure there are at least three unique row clusters
  clust_one <- col_clustering
  clust_two <- row_clustering
  # can I remove this one??
  if (n_clusts_row ==1 ) {
    clust_two <- cbind(clust_two, rbinom(nrow(row_clustering), 1, 0.1))
    n_clusts_row <- ncol(clust_two)
  }
  while (check(clust_two)) {
    clust_two <- cbind(clust_two, rbinom(nrow(row_clustering), 1, 0.1))
    n_clusts_row <- ncol(clust_two)
  }
  rep <- ifelse(n_clusts_row == n_clusts, FALSE, TRUE)

  # calculate score for each cluster
  s_vals <- vector("list", length=n_clusts)
  for (k in 1:n_clusts){
    indices <- clust_two[, k] == 1
    # if row or col cluster empty, set score to 0
    if ((sum(indices) == 0) | (sum(clust_one[,k] == 1) == 0)) {
      bisil_score[k] <- bisil_score[k] + 0
    } else {
      # subset data using column cluster k
      new_data <- Xinput[, (clust_one[, k] == 1)]
      distances <- as.matrix(stats::dist(new_data, method))
      b_vec <- c()
      # calculate a values
      if(sum(indices)==1){
        # if only one element belongs to row cluster, set a=0
        a_vals <- 0
      }else{
        a_vals <- apply(distances[indices, indices], 
                        1, function(x) sum(x) / (length(x) - 1))
      }
      # calculate b values
      # indices for other clusters
      other <- (1:n_clusts_row)[-k]
      b_vals <- vector("list", length = (n_clusts_row - 1))
      t <- 1
      # consider every other cluster
      for (l in other) {
        oth_ind <- clust_two[, l] == 1
        if ((sum(oth_ind) == 0) | all(oth_ind == indices)) {
          # if the other cluster is empty,
          # or if the other cluster is the same as the current cluster, set b=Inf
          b_val <- rep(Inf, sum(indices))
        } else if ((sum(oth_ind) == 1) | (sum(indices) == 1)){
          # if either cluster has only one element, need to calculate mean differently 
          b_val <- mean(distances[indices, oth_ind])
        } else {
          b_val <- rowMeans(distances[indices, oth_ind])
        }
        b_vec <- c(b_vec, mean(b_val))
        b_vals[[t]] <- b_val
        t <- t + 1
      }
      b_vals <- b_vals[[which.min(b_vec)]]
      # calculate bisilhouette score
      if (all(b_vals == Inf) | (all(b_vals == 0) & all(a_vals == 0))) {
        # if all b values are Inf, set score to 0
        # this corresponds to all other clusters being empty, or the same as the current cluster
        # if all a and b values are 0, set score to 0
        # this corresponds to the current cluster being empty
        s_vals[[k]] <- 0
        bisil_score[k] <- bisil_score[k] + 0
      } else {
        s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
        s_vals[[k]] <- s_con
        bisil_score[k] <- bisil_score[k] + mean(s_con)
      }
    }
  }
  # calculate overall score
  if (sum(bisil_score) == 0){
    bisil <- 0
  }else{
    bisil <- ifelse(sum(bisil_score != 0) == 1, sum(bisil_score),
                sum(bisil_score) / (sum(bisil_score != 0)) -
                2 * sd(bisil_score[bisil_score != 0]))
  }  
  #return relationships and mean of bisil_score across views
  return(list("bisil" = bisil, "vals" = s_vals, "repeat" = rep))
}



BisilhouetteScore <- function(Xinput, row_clustering, col_clustering, method="euclidean", seed=TRUE, n_reps=10){
   # Calculate the bisilhouette score  
    #
    # Args:
    #  Xinput: data matrix, shape (N, p).
    #  row_clustering: binary matrix indicating row clustering, shape(N, k).
    #  col_clustering: binary matrix indicating column clustering, shape(p, k).
    #  method: distance metric to use, str. Default is "euclidean".
    #
    # Returns:
    #  bisil: bisilhouette score, float.
 
  # Error handling
  if (nrow(Xinput) != nrow(row_clustering)) {
      stop("Number of rows in data and row clustering do not match.")
  }
  if (ncol(Xinput) != nrow(col_clustering)) {
      stop("Number of columns in data and column clustering do not match.")
  }
  if (ncol(row_clustering) != ncol(col_clustering)) {
      stop("Number of row and column clusters do not match.")
  }

  # set seed 
  if (!seed) {
    set.seed(seed)
  }

  # normalise input data if needed
  if (any(colSums(Xinput) != 1)) {
    Xinput <- single_alt_l1_normalisation(Xinput)$newMatrix
  }

  
  
  # initial results
  results <- BisilScoreInner(Xinput, row_clustering, col_clustering, method)
  bisil <- results$bisil
  vals <- results$vals
  # repeat if necessary 
  if(results$rep){
    for(i in 1:10){
        res_rep <- BisilScoreInner(Xinput, row_clustering, col_clustering, method)
        bisil <- bisil + res_rep$bisil
        vals <- lapply(1:length(vals), function(k) vals[[k]] + res_rep$vals[[k]])
    }
    bisil <- bisil / 10
    vals <- lapply(vals, function(x) x / 10)
  }
  return(list("bisil" = bisil, "vals" = vals))
}