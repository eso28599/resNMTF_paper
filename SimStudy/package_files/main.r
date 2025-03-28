# Libraries
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
library(bisilhouette)


# Application of ResNMTF to data!
#'
#' @param
rest_multi_nmtf_inner <- function(data, init_f = NULL, init_s = NULL,
                                  init_g = NULL, k_vec = NULL,
                                  phi = NULL, xi = NULL, psi = NULL,
                                  n_iters = NULL,
                                  repeats = 5,
                                  distance = "euclidean",
                                  no_clusts = FALSE) {
  #' Run Restrictive-Multi NMTF, following the above algorithm
  #'
  #' 1. Normalise X^v s.t ||X||_1 = 1
  #' 2. Initialise F, S, G
  #' 3. Repeat for each view u
  #'    a. Fix ALL, update F^u
  #'    c. Fix ALL, update G^u
  #'    d. Fix ALL, update S^u
  #' until convergence
  n_v <- length(data)

  # initialise F, S and G based on svd decomposition if not given
  if (is.null(init_f) | is.null(init_g) | is.null(init_s)) {
    inits <- init_mats(data, k_vec)
    current_f <- inits$init_f
    current_s <- inits$init_s
    current_g <- inits$init_g
    currentlam <- inits$init_lambda
    currentmu <- inits$init_mu
  } else {
    # Take init_f, init_s, init_g as the initialised latent representations
    current_f <- init_f
    current_s <- init_s
    current_g <- init_g
    currentlam <- lapply(current_f, colSums)
    currentmu <- lapply(current_g, colSums)
  }
  # Initialising additional parameters
  x_hat <- vector("list", length = n_v)
  # Update until convergence, or for n_iters times
  if (is.null(n_iters)) {
    total_err <- c()
    # Run while-loop until convergence
    err_diff <- 1
    err_temp <- 0
    while ((err_diff > 1.0e-6)) {
      err <- numeric(length = n_v)
      new_parameters <- update_matrices(
        X = data,
        init_f = current_f,
        init_s = current_s,
        init_g = current_g,
        lambda = currentlam,
        mu = currentmu,
        phi = phi,
        xi = xi,
        psi = psi
      )
      current_f <- new_parameters$output_f
      current_s <- new_parameters$output_s
      current_g <- new_parameters$output_g
      currentlam <- new_parameters$lamoutput
      currentmu <- new_parameters$muoutput
      for (v in 1:n_v) {
        x_hat[[v]] <- current_f[[v]] %*% current_s[[v]] %*% t(current_g[[v]])
        err[v] <- sum((data[[v]] - x_hat[[v]])**2) / sum((data[[v]])**2)
      }
      mean_err <- mean(err)
      total_err <- c(total_err, mean_err)
      err_diff <- abs(mean_err - err_temp)
      err_temp <- tail(total_err, n = 1)
    }
  } else {
    total_err <- numeric(length = n_iters)
    for (t in 1:n_iters) {
      err <- numeric(length = length(current_f))
      new_parameters <- update_matrices(
        X = data,
        init_f = current_f,
        init_s = current_s,
        init_g = current_g,
        lambda = currentlam,
        mu = currentmu,
        phi = phi,
        xi = xi,
        psi = psi
      )
      current_f <- new_parameters$output_f
      current_s <- new_parameters$output_s
      current_g <- new_parameters$output_g
      currentlam <- new_parameters$lamoutput
      currentmu <- new_parameters$muoutput
      for (v in 1:n_v) {
        x_hat[[v]] <- current_f[[v]] %*% current_s[[v]] %*% t(current_g[[v]])
        err[v] <- sum((data[[v]] - x_hat[[v]])**2) / sum((data[[v]])**2)
      }
      total_err[t] <- mean(err)
    }
  }
  for (v in 1:n_v) {
    normal_f <- matrix_normalisation(current_f[[v]])
    current_f[[v]] <- normal_f$normalised_matrix
    normal_g <- matrix_normalisation(current_g[[v]])
    current_g[[v]] <- normal_g$normalised_matrix
    current_s[[v]] <- (normal_f$normaliser) %*% current_s[[v]] %*%
      normal_g$normaliser
  }
  # if only need to obtain factorisation, return values now
  if (no_clusts) {
    return(list(
      "output_f" = current_f, "output_s" = current_s,
      "output_g" = current_g
    ))
  }
  # find clustering results and bisilhouette score
  clusters <- obtain_biclusters(
    data, current_f,
    current_g, current_s, repeats, distance
  )
  if (is.null(n_iters)) {
    error <- mean(tail(total_err, n = 10))
  } else {
    error <- tail(total_err, n = 1)
  }
  return(list(
    "output_f" = current_f, "output_s" = current_s,
    "output_g" = current_g, "Error" = error,
    "All_Error" = total_err, "Sil_score" = clusters$bisil,
    "row_clusters" = clusters$row_clustering,
    "col_clusters" = clusters$col_clustering,
    "lambda" = currentlam,
    "mu" = currentmu
  ))
}

#' @param k_max integer, default is 6, must be greater than 2, largest value of k to be considered initially,
#' @param k_min integer, default is 3, must be greater than 1, smallest value of k to be considered initially,
#' @param distance string, default is "euclidean", distance metric to use for clustering,
#' @param repeats integer, default is 5, minimum value of 2, number of repeats to use for clustering,
#' @param no_clusts boolean, default is FALSE, whether to return only the factorisation or not,
#' @param sample_rate numeric, default is 0.9, proportion of data to sample for stability analysis,
#' @param n_stability integer, default is 5, number of times to repeat stability analysis,
#' @param stability boolean, default is TRUE, whether to perform stability analysis or not,
#' @param stab_thres numeric, default is 0.4, threshold for stability analysis,
#' @param stab_test boolean, default is FALSE, whether to perform stability test or not,
#' @param data list of matrices, data to be factorised,
#' @param init_f list of matrices, initialisation for F matrices,
#' @param init_s list of matrices, initialisation for S matrices,
#' @param init_g list of matrices, initialisation for G matrices,
#' @param k_vec integer, vector of integers, number of clusters to consider,
#' @param phi list of matrices, default is NULL, restriction matrices for F,
#' @param xi list of matrices, default is NULL, restriction matrices for S,
#' @param psi list of matrices, default is NULL, restriction matrices for G,
#' @param n_iters integer, default is NULL, number of iterations to run for,
#' @return list of results from ResNMTF
#' @export
#' @examples
#' data <- list(matrix(rnorm(100), nrow = 10), matrix(rnorm(100), nrow = 10))
#' restMultiNMTF_run(data = data, k_vec = c(3, 3), n_iters = 100)
#' restMultiNMTF_run(
#'   data = data, k_vec = c(3, 3), n_iters = 100,
#'   stability = FALSE
#' )
#' restMultiNMTF_run(data = data, k_vec = c(3, 3), n_iters = 100, no_clusts = TRUE)
#' restMultiNMTF_run(data = data, k_vec = c(3, 3), n_iters = 100, k_min = 3, k_max = 8)
#' restMultiNMTF_run(data = data, k_vec = c(3, 3), n_iters = 100, k_min = 3, k_max = 8, distance = "euclidean")
#' restMultiNMTF_run(data = data, k_vec = c(3, 3), n_iters = 100, k_min = 3, k_max = 8, distance = "euclidean", repeats = 5)
#' restMultiNMTF_run(data = data, k_vec = c(3, 3), n_iters = 100, k_min = 3, k_max = 8, distance = "euclidean", repeats = 5, no_clusts = FALSE)
restMultiNMTF_run <- function(data, init_f = NULL, init_s = NULL,
                              init_g = NULL, k_vec = NULL,
                              phi = NULL, xi = NULL, psi = NULL,
                              n_iters = NULL, k_min = 3, k_max = 8, distance = "euclidean", repeats = 5, no_clusts = FALSE,
                              sample_rate = 0.9, n_stability = 5, stability = TRUE, stab_thres = 0.4, stab_test = FALSE) {
  # initialise phi etc matrices as zeros if not specified
  # otherwise multiply by given parameter
  n_v <- length(data)
  data <- make_non_neg(data)
  if (!typeof(data[[1]]) == "double") {
    data <- lapply(data, function(x) as.matrix(x))
  }
  # Normalise data
  data <- lapply(data, function(x) matrix_normalisation(x)$normalised_matrix)
  # initialise restriction matrices if not specified
  # views with no restrictions require no input
  phi <- init_rest_mats(phi, n_v)
  psi <- init_rest_mats(psi, n_v)
  xi <- init_rest_mats(xi, n_v)
  # if number of clusters has been specified method can be applied straight away
  if ((!is.null(k_vec))) {
    results <- rest_multi_nmtf_inner(
      data, init_f, init_s, init_g,
      k_vec, phi, xi, psi, n_iters,
      repeats, distance, no_clusts
    )
    # if using the original data, we want to perform stability analysis
    # otherwise we want the results
    if (stability) {
      return(stability_check(
        data, init_s, results,
        k_vec, phi, xi, psi, n_iters,
        repeats, no_clusts, distance, sample_rate,
        n_stability, stab_thres
      ))
    } else {
      return(results)
    }
  }
  # define set of k_s to consider
  k_vec <- k_min:k_max
  k_vec <- rep(1, n_v)
  n_k <- length(k_vec)
  # initialise storage of results
  # apply method for each k to be considered
  # how many jobs you want the computer to run at the same time
  # if on windows operating system - do normal for loop
  if ((.Platform$OS.type == "windows") | (.Platform$OS.type == "unix")) {
    res_list <- vector("list", length = n_k)
    for (i in 1:n_k) {
      res_list[[i]] <- rest_multi_nmtf_inner(
        data, init_f, init_s, init_g,
        k_vec[i] * k_vec, phi, xi, psi, n_iters,
        repeats, distance, no_clusts
      )
    }
  } else {
    # Get the total number of cores
    numberOfCores <- detectCores()
    # Register all the cores
    registerDoParallel(min(numberOfCores, length(k_vec)))
    res_list <- foreach(i = 1:length(k_vec)) %dopar% {
      rest_multi_nmtf_inner(
        data, init_f, init_s, init_g,
        k_vec[i] * k_vec, phi, xi, psi, n_iters,
        repeats, distance, no_clusts
      )
    }
  }
  # extract scores
  err_list <- rep(0, length(k_vec))
  for (i in 1:length(k_vec)) {
    err_list[i] <- res_list[[i]][["Sil_score"]][1]
  }
  # find value of k of lowest error
  test <- k_vec[which.max(err_list)]
  max_i <- k_max
  # if best performing k is the largest k considered
  # apply method to k + 1 until this is no longer the case
  if (k_min != k_max) {
    while (test == max_i) {
      max_i <- max_i + 1
      k_vec <- c(k_vec, max_i)
      k <- max_i * k_vec
      new_l <- length(k_vec)
      res_list[[new_l]] <- rest_multi_nmtf_inner(
        data, init_f, init_s, init_g,
        k, phi, xi, psi, n_iters,
        repeats, distance, no_clusts
      )
      err_list <- c(err_list, res_list[[new_l]][["Sil_score"]][1])
      test <- k_vec[which.max(err_list)]
    }
  }
  print(err_list)
  k <- which.max(err_list)
  results <- res_list[[k]]
  if (stability) {
    return(stability_check(
      data, init_s, results,
      k_vec[k], phi, xi, psi, n_iters,
      repeats, no_clusts, distance,
      sample_rate, n_stability,
      stab_thres, stab_test
    ))
  } else {
    return(results)
  }
}
