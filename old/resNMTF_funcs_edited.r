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


## General manipulation functions 
check_neg <- function(X){
  #ensures a list of matrices are non-negative
  return(lapply(X, make_non_neg))
}

make_non_neg <- function(X_mat){
  #makes a matrix non-negative
  return(apply(X_mat, 2, function(x) x + abs(min(0, min(x)))))
}


star_prod <- function(vec,mat_list){ 
  # function which takes a list L and vector v as input 
  # returns sum(v_i*L[[i]])
  #' mat_list: list of k matrices
  #' vec: a vector of length k
  #' Output: a matrix of same dimensions as the entries in mat_list
  vec_mat <- 0 
  for(i in 1:length(vec)){
    if(vec[i]!=0){
      vec_mat <- vec_mat + vec[i] * mat_list[[i]]
    }
  }
  return(vec_mat)
}

single_alt_l1_normalisation <- function(Xmatrix){
  # normalises a matrix so that the l1 norm of each column is 1
  Q <- diag(colSums(Xmatrix))
  #solve Q
  newMatrix <- Xmatrix %*% solve(Q)
  return(list("Q" = Q, "newMatrix" = newMatrix))
}

##Initialisation functions
init_rest_mats <- function(mat, n_v){
  # inititialise phi/psi/xi matrices
  # Check if the input matrix is NULL
  if (is.null(mat)){
    # If NULL
    # initialize a new matrix of zeros with dimensions determined by n_v
    return(matrix(0, nrow = n_v, ncol = n_v))
  } else {
    # If not NULL, return the input matrix (validate the existing matrix)
    return(mat)
  }
}


init_mats <- function(X, KK, sigma_I = 0.05){ 
    #' X: list of input data
    # Initialisation of F, S, G lists
    n_views <- length(X)
    Finit <- vector("list", length = n_views)
    Sinit <- vector("list", length = n_views)
    Ginit <- vector("list", length = n_views)
    lambda_init <- vector("list", length = n_views)
    mu_init <- vector("list", length = n_views)
    #initialise based on svd
    for(i in 1:n_views) {
        K <- KK[i]
        vals <- 1:K
        ss <- svd(X[[i]])
        Finit[[i]] <- abs(ss$u[, vals])
        F_normal <- single_alt_l1_normalisation(Finit[[i]])
        Finit[[i]] <- F_normal$newMatrix
        Sinit[[i]] <- abs(diag(ss$d)[vals, vals])
        Sinit[[i]] <- Sinit[[i]] + abs(mvrnorm(n = K,
                   mu = rep(0, K), Sigma = sigma_I * diag(K))[vals, vals])
        Ginit[[i]] <- abs(ss$v[, vals])
        G_normal <- single_alt_l1_normalisation(Ginit[[i]])
        Ginit[[i]] <- G_normal$newMatrix
        Sinit[[i]]  <- (F_normal$Q) %*% Sinit[[i]] %*% G_normal$Q
        lambda_init[[i]] <- colSums(Finit[[i]])
        mu_init[[i]] <-  colSums(Ginit[[i]])
    }

return(list("Finit" = Finit, "Ginit" = Ginit, "Sinit" = Sinit,
             "lambda_init" = lambda_init, "mu_init" = mu_init))
}



## Udate functions
update_F <- function(Xinput, Finput, Sinput, Ginput, lambda_in, phi, k){
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
  denominator_matrix <- currentF %*% Sinput %*% t(Ginput) %*%  Ginput %*% t(Sinput)

  #calculate the column vector based on phi that is needed
  phi_vec <- (phi + t(phi))[, k]
  lambda_mat <- 0.5 * t(matrix(lambda_in, length(lambda_in), nrow(currentF)))
  if (sum(phi_vec) == 0) {
    outputF <- currentF * ((numerator_matrix) / (denominator_matrix + lambda_mat))
  }else {
    num_mat_prod <- star_prod(phi_vec, Finput)
    denom_mat_prod <- sum(phi_vec) * currentF
    outputF <- currentF * ((numerator_matrix + num_mat_prod) / (denominator_matrix + denom_mat_prod + lambda_mat))
  }
  return(abs(outputF))
}

update_G <- function(Xinput, Finput, Sinput, Ginput,mu_in, psi, k){
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
  denominator_matrix <- currentG  %*% t(Sinput) %*% t(Finput) %*% Finput %*% Sinput
  mu_mat <- 0.5 * t(matrix(mu_in, length(mu_in), nrow(currentG)))
  #check if they are all the same size

  #calculate the column vector based on phi that is needed
  if (sum(psi) == 0) {
    outputG <- currentG * (numerator_matrix / (denominator_matrix + mu_mat))
  }else {
    psi_vec <- (psi + t(psi))[, k]
    num_mat_prod <- star_prod(psi_vec, Ginput)
    denom_mat_prod <- sum(psi_vec) * currentG
    outputG <- currentG * ((numerator_matrix + num_mat_prod) / (denominator_matrix + denom_mat_prod + mu_mat))
  }
  return(abs(outputG))
}

update_S <- function(Xinput, Finput, Sinput, Ginput, xi, k){
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

  #calculate the column vector based on phi that is needed
  if (sum(xi) == 0){
    outputS <- currentS * (numerator_matrix  / denominator_matrix)
  }else {
    xi_vec <- (xi + t(xi))[, k]
    num_mat_prod <- star_prod(xi_vec, Sinput)
    denom_mat_prod <- sum(xi_vec) * currentS
    outputS <- currentS * ((numerator_matrix + num_mat_prod) / (denominator_matrix + denom_mat_prod))
  }
  return(abs(outputS))
}

update_lm <- function(vec, matrix){
    #update lagrange multipliers
    return(colSums(matrix)* vec)
}

update_matrices <- function(Xinput, Finput, Sinput, Ginput, lambda, mu, phi, xi, psi, nIter){
  #' performs one interation of updates
  #' X: Input matrix
  #' phi: weight on restrictions for F -> matrix of size (n_v x n_v)
  #' psi: weight on restrictions for S -> matrix of size (n_v x n_v)
  #' xi: weight on restrictions for G -> matrix of size (n_v x n_v)
  #' Finit: Inital F matrix
  #' Sinit: Inital S matrix
  #' Ginit: Inital G matrix
  #' Output: Foutput, Soutput, Goutput
  
  # Update view-by-view
  n_v <- length(Xinput)
  currentF <- Finput
  currentS <- Sinput
  currentG <- Ginput
  currentlam <- lambda
  currentmu <- mu
  for (v in 1:n_v){
    # Update F
    currentF[[v]] <- update_F(Xinput = Xinput[[v]],
                              Finput = currentF,
                              Sinput = currentS[[v]],
                              Ginput = currentG[[v]],
                              lambda_in = currentlam[[v]],
                              phi = phi,
                              k = v)
    # Update G
    currentG[[v]] <- update_G(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS[[v]],
                              Ginput = currentG,
                              mu_in = currentmu[[v]],
                              psi = psi, k = v)
    # Update S
    currentS[[v]] <- update_S(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS,
                              Ginput = currentG[[v]],
                              xi = xi, k = v)
    currentlam[[v]] <-  update_lm(currentlam[[v]], currentF[[v]])
    currentmu[[v]] <-   update_lm(currentmu[[v]], currentG[[v]])
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  lam_output <- currentlam
  mu_output <- currentmu
  return(list("Foutput" = Foutput, "Soutput" = Soutput, "Goutput" = Goutput,
                "lamoutput" = lam_output, "muoutput" = mu_output))
}

##functions to remove spurious biclusters
###functions for calculating JSD 
jsd_calc <- function(x1, x2){
    #calculate the Jensen Shannon divergence between 
    #vectors x1 and x2
    max_val <- max(x1,x2)
    d1  <- density(x1, from=0, to=max_val)
    d2  <- density(x2, from=0, to=max_val)
    d1$y[d1$x>max(x1)] <- 0
    d2$y[d2$x>max(x2)] <- 0
    return(suppressMessages(JSD(rbind(d1$y, d2$y), unit = "log2", est.prob ="empirical")))
}

get_thresholds <- function(Xinput, Foutput, repeats){
    # shuffle data and 
    # find the threshold for removal of spurious biclusters
    #repeats: minimum value of 2
    n_views <- length(Xinput)
    n_clusts <- dim(Foutput[[1]])[2]
    k_input <- n_clusts * rep(1, length = n_views)
    x_mess <- vector(mode = "list", length = repeats)
    for (n in 1:repeats) {
      data_messed <- vector(mode = "list", length = n_views)
      for (i in 1:n_views){
          #correct shuffling
          dims <- dim(Xinput[[i]])
          data_messed[[i]] <- matrix(sample(Xinput[[i]]),
                                  dims[1], dims[2])
          while(any(colSums(data_messed[[i]])==0)| any(rowSums(data_messed[[i]])==0)){
              data_messed[[i]] <- matrix(sample(Xinput[[i]]),
                                       dims[1], dims[2])
          }
      }
      results <- restMultiNMTF_run(Xinput = data_messed,
               KK = k_input,  no_clusts = TRUE, stability = FALSE)
      x_mess[[n]] <- results$Foutput
    }
    avg_score <- c()
    max_score <- c()
    data <- vector(mode = "list", length = n_views)
    dens_list <- vector(mode = "list", length = n_views)
    d <- 1
    for (i in 1:n_views){
      scores <- c()
        for (j in 1:max(repeats - 1,1)){
          data[[i]] <- cbind(data[[i]], x_mess[[j]][[i]])
          for (k in 1:n_clusts){
            #jth repeat, ith view, kth cluster
            x1 <- ((x_mess[[j]])[[i]])[, k]
            for (l in (j + 1):repeats){
              for (m in 1:n_clusts){
                x2 <- ((x_mess[[l]])[[i]])[, m]
                max_val <- max(x1,x2)
                d1  <- density(x1, from=0, to=max_val)
                d2  <- density(x2, from=0, to=max_val)
                d1$y[d1$x>max(x1)] <- 0
                d2$y[d2$x>max(x2)] <- 0
                dens_list[[d]] <- d1
                scores <- c(scores,
                   suppressMessages(JSD(rbind(d1$y, d2$y), unit = "log2", est.prob="empirical")))
              }
            }
          }
      }
      data[[i]] <- cbind(data[[i]], x_mess[[repeats]][[i]])
      avg_score <- c(avg_score, mean(scores))
      #max_score <- c(max_score, max(scores))
      dens <- density(scores)
      max_score <- c(max_score, dens$x[which.max(dens$y)])
    }
  return(list("avg_score" = avg_score,
             "max_score" = max_score, "scores" = scores, "data" = data, "dens" = dens_list))
}


check_biclusters <- function(Xinput, Foutput, repeats){
  # calculate JSD between returned F and noise
  # returns scores, avg score and max score
  n_views <- length(Xinput)
  n_clusts <- dim(Foutput[[1]])[2]
  # updated results
  scores <- matrix(0, nrow = n_views, ncol = n_clusts)
  thresholds <- get_thresholds(Xinput, Foutput, repeats)
  for (i in 1:n_views){
    x_noise <- thresholds$data[[i]]
    for (k in 1:n_clusts){
      x <- Foutput[[i]][, k]
      scores[i, k] <- mean(apply(x_noise, 2, 
           function(y) jsd_calc(x,y)))
    }
  }
  return(list("score" = scores,
          "avg_threshold" = thresholds$avg_score,
          "max_threshold" = thresholds$max_score))
}

##Functions to return clustering with spurious biclusters removed

bisil_score <- function(Xinput, row_clustering, col_clustering, final=FALSE){
    #simultaneously calculates bisilhouette score for 
    #clustering as well matching clusters correctly.
    n_views <- length(Xinput)
    n_clusts <- dim(row_clustering[[1]])[2]
    relations <- matrix(0, nrow = n_views, n_clusts)
    sil_score <- matrix(0, nrow = n_views, n_clusts)
    clust_one <- col_clustering
    clust_two <- row_clustering
    for (i in 1:n_views){
        s_mat <- matrix(0, nrow = n_clusts, n_clusts)
        a_mat <- matrix(0, nrow = n_clusts, n_clusts)
        b_mat <- matrix(0, nrow = n_clusts, n_clusts)
        for (k in 1:n_clusts){
          #select data from specific column clustering
          if (sum(clust_one[[i]][, k] == 1) == 0){
            s_mat[k, ] <- 0
          }else{
            new_data <- Xinput[[i]][, (clust_one[[i]][, k] == 1)]
            spear_dists <- as.matrix(dist(new_data))/ (dim(new_data)[2])
            for (j in 1:n_clusts){
              indices <- clust_two[[i]][, j] == 1
              if (sum(indices) == 0) {
                s_mat[k, j] <- 0
              }else{
                a_vals <- apply(spear_dists[indices, indices], 1, function(x) sum(x)/(length(x)-1))
                #other clusts
                other <- (1:n_clusts)[-j]
                b_vec <- c()
                b_vals <- vector("list", length = (n_clusts - 1))
                t <- 1
                for(l in other){
                    oth_ind <- clust_two[[i]][, l] == 1
                    if(sum(oth_ind)==0){
                      b_val <- Inf
                    }else if (all(oth_ind==indices)) {
                       b_val <- rep(Inf, length(oth_ind))
                    }else{
                      if((sum(oth_ind)==1)|(sum(indices)==1)){
                        b_val <- mean(spear_dists[indices, oth_ind])
                      }else{
                        b_val <- rowMeans(spear_dists[indices, oth_ind])
                      }
                    }
                    b_vec <- c(b_vec, mean(b_val))
                    b_vals[[t]] <- b_val
                    t <- t + 1
                }
                closest <- which.min(b_vec)
                b_vals <- b_vals[[closest]]
                s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
                s_mat[k, j] <- mean(s_con)
               }
            }
          }
        }
        if(final){
          #if final clustering, don't want to change order
          relations[i, ] <- apply(s_mat, 1, which.min)
          sil_score[i, ] <- diag(s_mat)
        }else{
          relations[i, ] <- apply(s_mat, 1, which.max)
          sil_score[i, ] <- apply(s_mat, 1, max)
        }
    }
  overall_score <- apply(sil_score, 1, function(x) ifelse(sum(x) == 0,
               0, ifelse(sum(x!=0) == 1, mean(x[x!= 0]),
                 mean(x[x != 0]) - 2 * sd(x[x != 0]))))
  overall_score <- apply(sil_score, 1, function(x) ifelse(sum(x) == 0,
               0,mean(x[x!= 0])))
  overall_score <- ifelse(sum(overall_score)==0,
           0, mean(overall_score[overall_score!=0]))
  #return relationships and mean of sil_score across views
  return(list("scores" = sil_score, "relations" = relations, "overall" = overall_score))
}

clustering_res_NMTF <- function(Xinput, Foutput,
               Goutput, Soutput, repeats, epsilon = 0.02){
  #takes F,S,G returns clustering and removes 
  #spurious bicluster
  n_views <- length(Foutput)
  row_clustering <- vector("list", length = n_views)
  col_clustering <- vector("list", length = n_views)
  #check clustering and remove if necessary
  biclusts <- check_biclusters(Xinput, Foutput, repeats)
  for (i in 1:n_views) {
    row_clustering[[i]] <- apply(Foutput[[i]],
           2, function(x) as.numeric(x > (1 / dim(Foutput[[i]])[1])))

    col_clustering[[i]] <- apply(Goutput[[i]],
           2, function(x) as.numeric(x > (1 / dim(Goutput[[i]])[1])))
  }
  sil_score <- bisil_score(Xinput,
       row_clustering, col_clustering)
  s_score <- sil_score$scores
  #update realtions and
  #set biclusters that aren't strong enough to 0
  #and if bicluster is empty set row and cols to 0
  for (i in 1:n_views){
     indices <- (((biclusts$score[i, ]) < biclusts$max_threshold[i])
                   | ((biclusts$score[i, ]) == 0))
    
     relations <- sil_score$relations[i, ]
     new_indices <- indices[relations] # i==0, col cluster i isn't a bicluster
     s_score[i, new_indices] <- 0
     row_clustering[[i]] <- row_clustering[[i]][, relations]
     row_clustering[[i]][, new_indices] <- 0
     col_clustering[[i]][, new_indices] <- 0
  }
  #caluclate bisil for each view
  sil <- apply(s_score, 1, function(x) ifelse(sum(x) == 0, 0,  mean(x[x!= 0])))
  #calculate overall bisil
  sil <- ifelse(sum(sil)==0, 0, mean(sil[sil!=0]))
  return(list("row_clustering" = row_clustering,
      "col_clustering" = col_clustering, "sil" = sil))
}


#Application of ResNMTF to data!
restMultiNMTF_main <- function(Xinput, Finput = NULL, Sinput = NULL,
         Ginput = NULL, KK = NULL,
          phi = NULL, xi = NULL, psi = NULL,
          nIter = NULL, 
         repeats = 5,no_clusts = FALSE){
  #' Run Restrictive-Multi NMTF, following the above algorithm
  #' 
  #' 1. Normalise X^v s.t ||X||_1 = 1
  #' 2. Initialise F, S, G
  #' 3. Repeat for each view u
  #'    a. Fix ALL, update F^u
  #'    c. Fix ALL, update G^u
  #'    d. Fix ALL, update S^u
  #' until convergence
  n_v <- length(Xinput)
  # check in correct form
  if (!typeof(Xinput[[1]]) == "double") {
    Xinput <- lapply(Xinput, function(x) as.matrix(x))
    }
  # Normalise Xinput
  Xinput <- lapply(Xinput, function(x) single_alt_l1_normalisation(x)$newMatrix)
  
  # initialise F, S and G based on svd decomposition if not given
  if (is.null(Finput) | is.null(Ginput) | is.null(Sinput) ) {
        inits <- init_mats(Xinput, KK)
        currentF <- inits$Finit
        currentS <- inits$Sinit
        currentG <- inits$Ginit
        currentlam <- inits$lambda_init
        currentmu <- inits$mu_init
  }else {
    # Take Finit, Sinit, Ginit as the initialised latent representations
    currentF <- Finput
    currentS <- Sinput
    currentG <- Ginput
    currentlam <- lapply(currentF, colSums)
    currentmu <- lapply(currentG, colSums)
  }
  # Initialising additional parameters
  Xhat <- vector("list", length = n_v)
  # Update until convergence, or for nIter times
  if (is.null(nIter)){
    total_err <- c()
    # Run while-loop until convergence
    err_diff <- 1
    err_temp <- 0
    while ((err_diff > 1.0e-6)) {
      err <- numeric(length = n_v)
      new_parameters <- update_matrices(X = Xinput,
                                           Finput = currentF,
                                           Sinput = currentS,
                                           Ginput = currentG,
                                           lambda = currentlam,
                                           mu = currentmu,
                                           phi = phi,
                                           xi = xi,
                                           psi = psi)
      currentF <- new_parameters$Foutput
      currentS <- new_parameters$Soutput
      currentG <- new_parameters$Goutput
      currentlam <- new_parameters$lamoutput
      currentmu <- new_parameters$muoutput
      for (v in 1:n_v){
        Xhat[[v]] <- currentF[[v]] %*% currentS[[v]] %*% t(currentG[[v]])
        err[v] <- sum((Xinput[[v]] - Xhat[[v]])**2)/ sum((Xinput[[v]])**2)
      }
      mean_err <- mean(err)
      total_err <- c(total_err, mean_err)
      err_diff <- abs(mean_err - err_temp)
      err_temp <- tail(total_err, n = 1)
    }
  } else {
    total_err <- numeric(length = nIter)
    for (t in 1:nIter){
      err <- numeric(length = length(currentF))
      new_parameters <- update_matrices(X = Xinput,
                                           Finput = currentF,
                                           Sinput = currentS,
                                           Ginput = currentG,
                                           lambda = currentlam,
                                           mu = currentmu,
                                           phi = phi,
                                           xi = xi,
                                           psi = psi)
      currentF <- new_parameters$Foutput
      currentS <- new_parameters$Soutput
      currentG <- new_parameters$Goutput
      currentlam <- new_parameters$lamoutput
      currentmu <- new_parameters$muoutput
      for (v in 1:n_v){
        Xhat[[v]] <- currentF[[v]] %*% currentS[[v]] %*% t(currentG[[v]])
        err[v] <- sum((Xinput[[v]] - Xhat[[v]])**2)/ sum((Xinput[[v]])**2)
      }
      total_err[t] <- mean(err)
    }
  }
  for(v in 1:n_v){
        F_normal <- single_alt_l1_normalisation(currentF[[v]])
        currentF[[v]] <- F_normal$newMatrix
        G_normal <- single_alt_l1_normalisation(currentG[[v]])
        currentG[[v]] <- G_normal$newMatrix
        currentS[[v]] <- (F_normal$Q) %*% currentS[[v]] %*% G_normal$Q
  }

  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  lam_out <- currentlam
  mu_out <- currentmu
  # if only need to obtain factorisation, return values now
  if (no_clusts) {
    return(list("Foutput" = Foutput, "Soutput" = Soutput,
              "Goutput" = Goutput))
  }
  # find clustering results and silhouette score
  clusters <- clustering_res_NMTF(Xinput, Foutput,
               Goutput, Soutput, repeats, phi)
  if (is.null(nIter)) {
    error <- mean(tail(total_err, n = 10))
  }else {
    error <- tail(total_err, n = 1)
  }
  return(list("Foutput" = Foutput, "Soutput" = Soutput,
              "Goutput" = Goutput, "Error" = error,
              "All_Error" = total_err, "Sil_score" = clusters$sil,
              "row_clusters" = clusters$row_clustering,
              "col_clusters" = clusters$col_clustering, 
              "lambda" = lam_out, 
              "mu" = mu_out))
}


##functions for stability selection

jaccard <- function(a, b) {
    #calculate jaccard between two vectors
    intersection <- length(intersect(a, b))
    union <-  length(a) + length(b) - intersection
    if (union == 0){
      return(0)
    }else{
      return(intersection / union)
    }
}

cart_prod <- function(a, b) {
  #returns cartesian product of two sets
  prod <- c()
  # check a or b are not empty sets
  if (length(a) == 0 || length(b) == 0) {
    return(NULL)
  }else{
    for (k in 1:length(a)){
      prod <- c(prod, paste(a[k], b))
    }
  return(prod)
  }
}
jaccard_results <- function(row_c, col_c, true_r, true_c, stability = FALSE){
  m <- ncol(row_c)
  n <- ncol(true_r)
  # if no biclusters detected but some are present
  # return 0
  # if no biclusters present but some are detected - score of 0
  m_0 <- sum(colSums(row_c) != 0) #no of clusters actually detected
  n_0 <- sum(colSums(true_r) != 0) #no of true clusters
  if ((m_0 == 0 && n_0 != 0) || (n_0 == 0 && m_0 != 0)) {
    if (stability) {
      return(0)
    }else{
      return(list("rec" = rep(0, 2), "rel" = rep(0, 2), "f_score" = rep(0, 2)))
    }
  }
  # if no biclusters present and none detected - score of 1
  if (m_0 == 0 && n_0 == 0) {
    if (stability) {
      return(1)
    }else{
      return(list("rec" = rep(1, 2), "rel" = rep(1, 2), "f_score" = rep(1, 2)))
    }
  }
  samps <- 1:nrow(row_c)
  feats <- 1:nrow(col_c)
  #initialise storage of jaccard index between pairs
  jac_mat <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m){
    r_i <- samps[row_c[, i] == 1]
    c_i <- feats[col_c[, i] == 1]
    m_i <- cart_prod(r_i, c_i)
    for (j in 1:n){
        tr_i <- samps[true_r[, j] == 1]
        tc_i <- feats[true_c[, j] == 1]
        m_j <- cart_prod(tr_i, tc_i)
        jac_mat[i, j] <- jaccard(m_i, m_j)
    }
  }
  if (stability) {
    return(apply(jac_mat, 2, max))
  }
  rel <- ifelse(sum(apply(jac_mat, 1, max) != 0) == 0, 0,
    sum(apply(jac_mat, 1, max)) / m_0)
  rec <- ifelse(sum(apply(jac_mat, 2, max) != 0) == 0, 0, 
    sum(apply(jac_mat, 2, max)) / n_0)
  f <- ifelse(rel * rec == 0, 0, 2 * rel * rec / (rel + rec))
  return(list("rec" = rep(rec, 2), "rel" = rep(rel, 2), "f_score" = rep(f, 2)))
}

test_cond <- function(data, attempt){
  if(attempt==1){
    return(TRUE)
  }
  return(any(unlist(lapply(data, function(x) any(colSums(x)==0) | any(rowSums(x)==0)))))
}

#stab_thres changed from 0.6 to 0.4
stability_check <- function(Xinput, results,
                     k, phi, xi, psi, nIter,
                     repeats, no_clusts, sample_rate = 0.9,
                     n_stability = 5, stab_thres = 0.4){
    #check whether stability check even needs to be done
    #no_clusts_detected
    n_c <- sum(as.numeric(lapply(results$row_clusters,
             function(x) sum(colSums(x)))))
    if (n_c == 0){
      print("No biclusters detected!")
      return(results)
    }
    n_views <- length(Xinput)
    # initialise storage of results
    jacc <- matrix(0, nrow = n_views, ncol = k[1])
    jacc_rand <- matrix(0, nrow = n_views, ncol = k[1])
    for (t in 1:n_stability){
      new_data <- vector(mode = "list", length = n_views)
      row_samples <- vector(mode = "list", length = n_views)
      col_samples <- vector(mode = "list", length = n_views)
      #turn this into a function to be used with lapply

      dim <- dim(Xinput[[1]])
      attempt <- 1
      while(test_cond(new_data, attempt)){
        if(attempt==20){
          print("Unable to perform stability analysis due to sparsity of data.")
          return(results)
        }
        row_samples[[1]] <- sample(dim[1], (dim[1] * sample_rate))
        col_samples[[1]] <- sample(dim[2], (dim[2] * sample_rate))
        new_data[[1]] <- Xinput[[1]][row_samples[[1]], col_samples[[1]]]
        if(any(colSums(new_data[[1]])==0) | any(rowSums(new_data[[1]])==0)){
              zeros_cols <- colSums(new_data[[1]])!=0
              zeros_rows <- rowSums(new_data[[1]])!=0
              row_samples[[1]] <- row_samples[[1]][zeros_rows]
              col_samples[[1]] <- col_samples[[1]][zeros_cols]
              new_data[[1]] <- Xinput[[1]][row_samples[[1]], col_samples[[1]]]
        }
        if(n_views>1){
          for(i in 2:n_views){
            dims <- dim(Xinput[[i]])
            if((dims[1])==dim[1]){
              row_samples[[i]] <- row_samples[[1]]
            }else{
              row_samples[[i]] <- sample(dims[1], (dims[1] * sample_rate))
            }
            if((dims[2])==dim[2]){
              col_samples[[i]] <- col_samples[[1]]
            }else{
              col_samples[[i]] <- sample(dims[2], (dims[2] * sample_rate))
            }
            new_data[[i]] <- Xinput[[i]][row_samples[[i]], col_samples[[i]]]
            if(any(colSums(new_data[[i]])==0) | any(rowSums(new_data[[i]])==0)){
              zeros_cols <- colSums(new_data[[i]])!=0
              zeros_rows <- rowSums(new_data[[i]])!=0
              if((dims[1])==dim[1]){
                for(p in 1:i){
                  row_samples[[p]] <- row_samples[[p]][zeros_rows]
                }
              }else{
                row_samples[[i]] <- row_samples[[i]][zeros_rows]
              }
              if((dims[2])==dim[2]){
                for(p in 1:i){
                  col_samples[[p]] <- col_samples[[p]][zeros_cols]
                }
              }else{
                col_samples[[i]] <- col_samples[[i]][zeros_cols]
              }
              new_data[[p]] <- Xinput[[p]][row_samples[[p]], col_samples[[p]]]
            }
        }
      }
      attempt <- attempt + 1
      }
      new_results <- restMultiNMTF_main(new_data, Finput = NULL, Sinput = NULL,
          Ginput = NULL, k,
          phi, xi, psi, nIter,repeats)
      #compare results
      #extract results
      for(i in 1:n_views){
        jacc[i, ] <- jacc[i, ] + jaccard_results(new_results$row_clusters[[i]],
                 new_results$col_clusters[[i]],
              results$row_clusters[[i]][row_samples[[i]], ],
               results$col_clusters[[i]][col_samples[[i]], ], TRUE)
        #jacc_rand[i, ] <- jacc_rand[i, ] + jaccard_results(apply(new_results$row_clusters[[i]], 2, sample),
         #        apply(new_results$col_clusters[[i]], 2, sample),
          #    results$row_clusters[[i]][row_samples[[i]], ],
           #    results$col_clusters[[i]][col_samples[[i]],], TRUE)
        jacc_rand[i, ] <- jacc_rand[i, ] + jaccard_results(new_results$row_clusters[[i]],
                 apply(new_results$col_clusters[[i]], 2, sample),
              results$row_clusters[[i]][row_samples[[i]], ],
               results$col_clusters[[i]][col_samples[[i]],], TRUE)
      }
   
    }
    jacc <- jacc / n_stability
    jacc_rand <- jacc_rand / n_stability
    for (i in 1:n_views){
      #set clusters not deemed stable to have 0 members
      results$row_clusters[[i]][, jacc[i, ] <  stab_thres] <- 0
      results$col_clusters[[i]][, jacc[i, ] <  stab_thres] <- 0
      results$row_clusters[[i]][, jacc[i, ] <  jacc_rand[i, ]] <- 0
      results$col_clusters[[i]][, jacc[i, ] <  jacc_rand[i, ]] <- 0
    }
    results$Sil_score <- bisil_score(Xinput,
                  results$row_clusters, results$col_clusters, TRUE)$overall
    return(results)
}


restMultiNMTF_run <- function(Xinput, Finput=NULL, Sinput=NULL, 
            Ginput=NULL, KK=NULL, phi=NULL, xi=NULL, psi=NULL, 
            nIter=NULL, k_min=3, k_max =6, repeats = 5, no_clusts = FALSE, 
             sample_rate = 0.8, n_stability = 5, stability = TRUE, stab_thres = 0.6){

  #' @param k_max integer, default is 6, must be greater than 2, largest value of k to be considered initially,
  # initialise phi etc matrices as zeros if not specified
  # otherwise multiply by given parameter
  n_v <- length(Xinput)
  Xinput <- check_neg(Xinput)
  # initialise restriction matrices if not specified 
  # views with no restrictions require no input
  phi <- init_rest_mats(phi, n_v)
  psi <- init_rest_mats(psi, n_v)
  xi <- init_rest_mats(xi, n_v)
  # if number of clusters has been specified method can be applied straight away
  if ((!is.null(KK))) {
    results <- restMultiNMTF_main(Xinput, Finput, Sinput, Ginput,
                     KK, phi, xi, psi, nIter,
                      repeats, no_clusts)
    # if using the original data, we want to perform stability analysis 
    # otherwise we want the results
    if (stability) {
      return(stability_check(Xinput, results,
                     KK, phi, xi, psi, nIter,
                    repeats, no_clusts, sample_rate,
                     n_stability, stab_thres))
    }else {
      return(results)
    }
  }
  # define set of k_s to consider
  KK <- k_min:k_max
  k_vec <- rep(1, n_v)
  n_k <- length(KK)
  #initialise storage of results
  # apply method for each k to be considered
  # how many jobs you want the computer to run at the same time
  #if on windows operating system - do normal for loop
  if (.Platform$OS.type == "windows"){
    res_list <- vector("list", length = n_k)
    for (i in 1:n_k){
      res_list[[i]] <- restMultiNMTF_main(Xinput, Finput, Sinput, Ginput,
                     KK[i] * k_vec, phi, xi, psi, nIter,
                     repeats, no_clusts)
    }
  }else{
    # Get the total number of cores
    numberOfCores <- detectCores()
    # Register all the cores
    registerDoParallel(min(numberOfCores, length(KK)))
    res_list <- foreach(i = 1:length(KK)) %dopar% {
    restMultiNMTF_main(Xinput, Finput, Sinput, Ginput,
                     KK[i] * k_vec, phi, xi, psi, nIter,
                     repeats, no_clusts)
  }
  }

  #extract scores
  err_list <- rep(0, length(KK))
  for (i in 1:length(KK)){
    err_list[i] <- res_list[[i]][["Sil_score"]][1]
  }
  #print(err_list)
  #find value of k of lowest error
  test <- KK[which.max(err_list)]
  max_i <- k_max
  # if best performing k is the largest k considered
  # apply method to k + 1 until this is no longer the case
  if(k_min != k_max){
    while (test == max_i) {
    max_i <- max_i + 1
    KK <- c(KK, max_i)
    k <- max_i * k_vec
    new_l <- length(KK)
    res_list[[new_l]] <- restMultiNMTF_main(Xinput, Finput, Sinput, Ginput,
                     k, phi, xi, psi, nIter,
                     repeats, no_clusts)
    err_list <- c(err_list, res_list[[new_l]][["Sil_score"]][1])
    test <- KK[which.max(err_list)]
  }
  }
  k <- which.max(err_list)
  results <- res_list[[k]]
  k_vec <- rep(KK[k], length = n_v)
  if(stability){
      return(stability_check(Xinput, results,
                     k_vec, phi, xi, psi, nIter,
                    repeats, no_clusts, sample_rate, n_stability, stab_thres))
  }else{
        return(results)
  }
}
