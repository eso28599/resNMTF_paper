# utils
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
## General manipulation functions
make_non_neg <- function(X) {
  # makes a list of matrices non-negative
  #' X: list of matrices
  return(lapply(X, make_non_neg_inner))
}

make_non_neg_inner <- function(X_mat) {
  # makes a matrix non-negative
  return(apply(X_mat, 2, function(x) x + abs(min(0, min(x)))))
}

star_prod <- function(vec, mat_list) {
  # function which takes a list L and vector v as input
  # returns sum(v_i*L[[i]])
  #' mat_list: list of k matrices
  #' vec: a vector of length k
  #' Output: a matrix of same dimensions as the entries in mat_list
  vec_mat <- 0
  for (i in 1:length(vec)) {
    if (vec[i] != 0) {
      vec_mat <- vec_mat + vec[i] * mat_list[[i]]
    }
  }
  return(vec_mat)
}

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
