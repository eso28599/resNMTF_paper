## General manipulation functions
make_non_neg <- function(X) {
  # makes a list of matrices non-negative
  #' X: list of matrices
  return(lapply(X, make_non_neg_inner))
}

make_non_neg_inner <- function(matrix) {
  # makes a matrix non-negative
  return(apply(matrix, 2, function(x) x + abs(min(0, min(x)))))
}

star_prod <- function(vec, mat_list) {
  # function which takes a list L and vector v as input
  # returns sum(v_i*L[[i]])
  #' mat_list: list of k matrices
  #' vec: a vector of length k
  #' Output: a matrix of same dimensions as the entries in mat_list
  vec_mat <- 0
  for (i in seq_len(length(vec))) {
    if (vec[i] != 0) {
      vec_mat <- vec_mat + vec[i] * mat_list[[i]]
    }
  }
  return(vec_mat)
}

matrix_normalisation <- function(matrix) {
  # normalises a matrix so that the l1 norm of each column is 1
  normaliser <- diag(colSums(matrix))
  # solve Q
  normalised_matrix <- matrix %*% solve(normaliser)
  return(list(
    "normaliser" = normaliser,
    "normalised_matrix" = normalised_matrix
  ))
}

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
