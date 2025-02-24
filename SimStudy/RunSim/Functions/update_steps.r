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
## Initialisation functions
init_rest_mats <- function(mat, n_v) {
  # inititialise phi/psi/xi matrices
  # Check if the input matrix is NULL
  if (is.null(mat)) {
    # If NULL
    # initialize a new matrix of zeros with dimensions determined by n_v
    return(matrix(0, nrow = n_v, ncol = n_v))
  } else {
    # If not NULL, return the input matrix (validate the existing matrix)
    return(mat)
  }
}

init_mats <- function(X, KK, sigma_I = 0.05) {
  #' X: list of input data
  # Initialisation of F, S, G lists
  n_views <- length(X)
  init_f <- vector("list", length = n_views)
  init_s <- vector("list", length = n_views)
  init_g <- vector("list", length = n_views)
  init_lambda <- vector("list", length = n_views)
  init_mu <- vector("list", length = n_views)
  # initialise based on svd
  for (i in 1:n_views) {
    K <- KK[i]
    vals <- 1:K
    ss <- svd(X[[i]])
    init_f[[i]] <- abs(ss$u[, vals])
    normal_f <- single_alt_l1_normalisation(init_f[[i]])
    init_f[[i]] <- normal_f$newMatrix
    init_s[[i]] <- abs(diag(ss$d)[vals, vals])
    init_s[[i]] <- init_s[[i]] + abs(mvrnorm(
      n = K,
      mu = rep(0, K), Sigma = sigma_I * diag(K)
    )[vals, vals])
    init_g[[i]] <- abs(ss$v[, vals])
    G_normal <- single_alt_l1_normalisation(init_g[[i]])
    init_g[[i]] <- G_normal$newMatrix
    init_s[[i]] <- (normal_f$Q) %*% init_s[[i]] %*% G_normal$Q
    init_lambda[[i]] <- colSums(init_f[[i]])
    init_mu[[i]] <- colSums(init_g[[i]])
  }

  return(list(
    "init_f" = init_f, "init_g" = init_g, "init_s" = init_s,
    "init_lambda" = init_lambda, "init_mu" = init_mu
  ))
}

init_mats_random <- function(X, K) {
  #' X: list of input data
  # Initialisation of F, S, G lists
  n_views <- length(X)
  init_f <- vector("list", length = n_views)
  init_s <- vector("list", length = n_views)
  init_g <- vector("list", length = n_views)
  init_lambda <- vector("list", length = n_views)
  init_mu <- vector("list", length = n_views)
  # initialise based on svd
  for (i in 1:n_views) {
    # ss <- svd(X[[i]])
    init_f[[i]] <- abs(mvrnorm(
      n = nrow(X[[i]]),
      mu = rep(0, K), Sigma = diag(K)
    ))
    normal_f <- single_alt_l1_normalisation(init_f[[i]])
    init_f[[i]] <- normal_f$newMatrix
    init_s[[i]] <- abs(mvrnorm(
      n = K,
      mu = rep(0, K), Sigma = diag(K)
    ))
    init_g[[i]] <- abs(mvrnorm(
      n = ncol(X[[i]]),
      mu = rep(0, K), Sigma = diag(K)
    ))
    G_normal <- single_alt_l1_normalisation(init_g[[i]])
    init_g[[i]] <- G_normal$newMatrix
    # init_s[[i]]  <- (normal_f$Q) %*% init_s[[i]] %*% G_normal$Q
    init_lambda[[i]] <- colSums(init_f[[i]])
    init_mu[[i]] <- colSums(init_g[[i]])
  }

  return(list(
    "init_f" = init_f, "init_g" = init_g, "init_s" = init_s,
    "init_lambda" = init_lambda, "init_mu" = init_mu
  ))
}

## Udate functions
update_F <- function(Xinput, Finput, Sinput, Ginput, lambda_in, phi, k) {
  #' X: Input matrix
  #' F: row-clustering -- Entire list as input of length n_v
  #' S: connection between row-clustering and column-clustering
  #' G: column-clustering
  #' phi: weight on restrictions for F -> matrix of size (n_v x n_v) - a sparse matrix
  #' k: which view to update
  #' Output: An update for Finput[[k]]

  # Find numerator
  currentF <- Finput[[k]]
  numerator_matrix <- Xinput %*% Ginput %*% t(Sinput)
  denominator_matrix <- currentF %*% Sinput %*% t(Ginput) %*% Ginput %*% t(Sinput)

  # calculate the column vector based on phi that is needed
  phi_vec <- (phi + t(phi))[, k]
  lambda_mat <- 0.5 * t(matrix(lambda_in, length(lambda_in), nrow(currentF)))
  if (sum(phi_vec) == 0) {
    outputF <- currentF * ((numerator_matrix) / (denominator_matrix + lambda_mat))
  } else {
    num_mat_prod <- star_prod(phi_vec, Finput)
    denom_mat_prod <- sum(phi_vec) * currentF
    outputF <- currentF * ((numerator_matrix + num_mat_prod) / (denominator_matrix + denom_mat_prod + lambda_mat))
  }
  return(abs(outputF))
}

update_G <- function(Xinput, Finput, Sinput, Ginput, mu_in, psi, k) {
  #' X: Input matrix
  #' F: row-clustering
  #' S: connection between row-clustering and column-clustering
  #' G: column-clustering -- Entire list as input of length n_v
  #' xi: weight on restrictions for G -> matrix of size (n_v x n_v) - a sparse matrix
  #' k: which view to update
  #' Output: An update for Ginput[[k]]

  # Find numerator
  currentG <- Ginput[[k]]
  numerator_matrix <- t(Xinput) %*% Finput %*% Sinput
  denominator_matrix <- currentG %*% t(Sinput) %*% t(Finput) %*% Finput %*% Sinput
  mu_mat <- 0.5 * t(matrix(mu_in, length(mu_in), nrow(currentG)))
  # check if they are all the same size

  # calculate the column vector based on phi that is needed
  if (sum(psi) == 0) {
    outputG <- currentG * (numerator_matrix / (denominator_matrix + mu_mat))
  } else {
    psi_vec <- (psi + t(psi))[, k]
    num_mat_prod <- star_prod(psi_vec, Ginput)
    denom_mat_prod <- sum(psi_vec) * currentG
    outputG <- currentG * ((numerator_matrix + num_mat_prod) / (denominator_matrix + denom_mat_prod + mu_mat))
  }
  return(abs(outputG))
}

update_S <- function(Xinput, Finput, Sinput, Ginput, xi, k) {
  #' X: Input matrix
  #' F: row-clustering
  #' S: connection between row-clustering and column-clustering -- Entire list as input of length n_v
  #' G: column-clustering
  #' xi: weight on restrictions for G -> matrix of size (n_v x n_v) - a sparse matrix
  #' k: which view to update
  #' Output: An update for Sinput[[k]]

  # Find numerator
  currentS <- Sinput[[k]]
  numerator_matrix <- t(Finput) %*% Xinput %*% Ginput
  denominator_matrix <- t(Finput) %*% Finput %*% currentS %*% t(Ginput) %*% Ginput

  # calculate the column vector based on phi that is needed
  if (sum(xi) == 0) {
    outputS <- currentS * (numerator_matrix / denominator_matrix)
  } else {
    xi_vec <- (xi + t(xi))[, k]
    num_mat_prod <- star_prod(xi_vec, Sinput)
    denom_mat_prod <- sum(xi_vec) * currentS
    outputS <- currentS * ((numerator_matrix + num_mat_prod) / (denominator_matrix + denom_mat_prod))
  }
  return(abs(outputS))
}

update_lm <- function(vec, matrix) {
  # update lagrange multipliers
  return(colSums(matrix) * vec)
}

update_matrices <- function(Xinput, Finput, Sinput, Ginput, lambda, mu, phi, xi, psi, nIter) {
  #' performs one interation of updates
  #' X: Input matrix
  #' phi: weight on restrictions for F -> matrix of size (n_v x n_v)
  #' psi: weight on restrictions for S -> matrix of size (n_v x n_v)
  #' xi: weight on restrictions for G -> matrix of size (n_v x n_v)
  #' init_f: Inital F matrix
  #' init_s: Inital S matrix
  #' init_g: Inital G matrix
  #' Output: Foutput, Soutput, Goutput

  # Update view-by-view
  n_v <- length(Xinput)
  currentF <- Finput
  currentS <- Sinput
  currentG <- Ginput
  currentlam <- lambda
  currentmu <- mu
  for (v in 1:n_v) {
    # Update F
    currentF[[v]] <- update_F(
      Xinput = Xinput[[v]],
      Finput = currentF,
      Sinput = currentS[[v]],
      Ginput = currentG[[v]],
      lambda_in = currentlam[[v]],
      phi = phi,
      k = v
    )
    # Update G
    currentG[[v]] <- update_G(
      Xinput = Xinput[[v]],
      Finput = currentF[[v]],
      Sinput = currentS[[v]],
      Ginput = currentG,
      mu_in = currentmu[[v]],
      psi = psi, k = v
    )
    # Update S
    currentS[[v]] <- update_S(
      Xinput = Xinput[[v]],
      Finput = currentF[[v]],
      Sinput = currentS,
      Ginput = currentG[[v]],
      xi = xi, k = v
    )
    currentlam[[v]] <- update_lm(currentlam[[v]], currentF[[v]])
    currentmu[[v]] <- update_lm(currentmu[[v]], currentG[[v]])
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  lam_output <- currentlam
  mu_output <- currentmu
  return(list(
    "Foutput" = Foutput, "Soutput" = Soutput, "Goutput" = Goutput,
    "lamoutput" = lam_output, "muoutput" = mu_output
  ))
}
