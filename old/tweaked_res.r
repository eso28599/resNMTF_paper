source("Functions/utils.r")
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
#bisilhouette functions

sil_score_inner <- function(Xinput, row_clustering, col_clustering, method="euclidean"){
    #simultaneously calculates silhouette score for a 
    #clustering as well matching clusters correctly.
    if(any(colSums(Xinput)!=1)){
      Xinput <- single_alt_l1_normalisation(Xinput)$newMatrix
    }
    n_clusts <- ncol(row_clustering)
    n_clusts_row <- n_clusts
    sil_score <- rep(0, length = n_clusts)
    clust_one <- col_clustering
    clust_two <- row_clustering
    if(n_clusts_row==1){
      clust_two <- cbind(clust_two, rbinom(nrow(row_clustering), 1, 0.1))
      n_clusts_row <- ncol(clust_two)
    }
    while(check(clust_two)){
      clust_two <- cbind(clust_two, rbinom(nrow(row_clustering), 1, 0.1))
      n_clusts_row <- ncol(clust_two)
    }
    # define whether repeat needs to happen
    rep <- ifelse(n_clusts_row == n_clusts, FALSE, TRUE)
    s_vals <- vector("list", length=n_clusts)
    for (k in 1:n_clusts){
        indices <- clust_two[, k] == 1
       # if row or col cluster empty, set to 0
        if ((sum(indices) == 0) | (sum(clust_one[,k] == 1) == 0)){
          sil_score[k] <- 0
        }else {
          #select data from specific column clustering
        new_data <- Xinput[, (clust_one[, k] == 1)]
        spear_dists <- as.matrix(stats::dist(new_data, method))
        b_vec <- c()
        #if only one element in row clust
        if(sum(indices)==1){
          a_vals <- 0
        }else{
          a_vals <- apply(spear_dists[indices, indices], 
                         1, function(x) sum(x)/(length(x)-1))
        }
        # calculate b values
        other <- (1:n_clusts_row)[-k]
        b_vals <- vector("list", length = (n_clusts_row - 1))
        t <- 1
        for(l in other){
            oth_ind <- clust_two[, l] == 1
            if(sum(oth_ind)==0){
              #if other is empty
              b_val <- rep(Inf, sum(indices))
            }else if(all(oth_ind == indices)){
              #if other is equal
              b_val <- rep(Inf, sum(indices))
            }else if((sum(oth_ind) == 1) | (sum(indices) == 1)){
              #if either has only one element
              b_val <- mean(spear_dists[indices, oth_ind])
            }else{
              b_val <- rowMeans(spear_dists[indices, oth_ind])
            }
            b_vec <- c(b_vec, mean(b_val))
            b_vals[[t]] <- b_val
            t <- t + 1
            }
          closest <- which.min(b_vec)
          b_vals <- b_vals[[closest]]
        if(all(b_vals == Inf)){
          s_vals[[k]] <- 0
          sil_score[k] <- 0
        }else if(all(b_vals == 0) & all(a_vals == 0)){
          sil_score[k] <- 0
        }else{
          s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
          s_vals[[k]] <- s_con
          sil_score[k] <- mean(s_con)
        }
    }
    }
  if (sum(sil_score) == 0){
    sil <- 0
  }else{
    sil <- ifelse(sum(sil_score != 0) == 1, sum(sil_score),
                 sum(sil_score) / (sum(sil_score != 0)) -
                 2 * sd(sil_score[sil_score != 0]))
    sil <- ifelse(sum(sil_score != 0) == 1, sum(sil_score),
                 sum(sil_score) / (sum(sil_score != 0)))
  }
  #return relationships and mean of sil_score across views
  return(list("sil" = sil, "vals" = s_vals, "repeat" = rep))
}



sil_score <- function(Xinput, row_clustering, col_clustering, method="euclidean", seed=TRUE){
  #for visualisation - set seed
  if(!seed){
    set.seed(seed)
  }
  #initial results
  results <- sil_score_inner(Xinput, row_clustering, col_clustering, method)
  sil <- results$sil
  vals <- results$vals
  #repeat if necessary 
  if(results$rep){
    for(i in 1:10){
        res_rep <- sil_score_inner(Xinput, row_clustering, col_clustering, method)
        sil <- sil + res_rep$sil
        vals <- lapply(1:length(vals), function(k) vals[[k]] + res_rep$vals[[k]])
    }
    sil <- sil/10
    vals <- lapply(vals, function(x) x/10)
  }
  return(list("sil" = sil, "vals" = vals))
}

clustering_res_NMTF <- function(Xinput, Foutput,
               Goutput, Soutput, repeats, epsilon=0.2){
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
  sil <- c()
  #update realtions and
  #set biclusters that aren't strong enough to 0
  #and if bicluster is empty set row and cols to 0
  for (i in 1:n_views){
     indices <- (((biclusts$score[i, ]) < biclusts$max_threshold[i])
                   | ((biclusts$score[i, ]) == 0))
     relations <- apply(Soutput[[i]], 1, which.max)
     new_indices <- indices[relations] # i==0, col cluster i isn't a bicluster
     row_clustering[[i]] <- row_clustering[[i]][, relations]
     row_clustering[[i]][, new_indices] <- 0
     col_clustering[[i]][, new_indices] <- 0
     sil <- c(sil,
       sil_score(Xinput[[i]], row_clustering[[i]], col_clustering[[i]])$sil)
  }
  #calculate overall bisil
  sil <- ifelse(sum(sil)==0, 0, mean(sil[sil!=0]))
  return(list("row_clustering" = row_clustering,
      "col_clustering" = col_clustering, "sil" = sil))
}

#only not calculating bisil
stability_check <- function(Xinput, results,
                     k, phi, xi, psi, nIter,
                     repeats, no_clusts, sample_rate = 0.9,
                     n_stability = 5, stab_thres = 0.6){
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
    # results$Sil_score <- bisil_score(Xinput,
    #               results$row_clusters, results$col_clusters, TRUE)$overall
    return(results)
}




## recent bisil checking
sil_score_inner <- function(Xinput, row_clustering, col_clustering, method="euclidean"){
    #simultaneously calculates silhouette score for a 
    #clustering as well matching clusters correctly.
    # if(any(colSums(Xinput)!=1)){
    #   Xinput <- single_alt_l1_normalisation(Xinput)$newMatrix
    # }
    n_clusts <- ncol(row_clustering)
    n_clusts_row <- n_clusts
    sil_score <- rep(0, length = n_clusts)
    clust_one <- col_clustering
    clust_two <- row_clustering
    if(n_clusts_row==1){
      clust_two <- cbind(clust_two, rbinom(nrow(row_clustering), 1, 0.1))
      n_clusts_row <- ncol(clust_two)
    }
    while(check(clust_two)){
      clust_two <- cbind(clust_two, rbinom(nrow(row_clustering), 1, 0.1))
      n_clusts_row <- ncol(clust_two)
    }
    # define whether repeat needs to happen
    rep <- ifelse(n_clusts_row == n_clusts, FALSE, TRUE)
    s_vals <- vector("list", length=n_clusts)
    for (k in 1:n_clusts){
        indices <- clust_two[, k] == 1
       # if row or col cluster empty, set to 0
        if ((sum(indices) == 0) | (sum(clust_one[,k] == 1) == 0)){
          sil_score[k] <- 0
        }else {
          #select data from specific column clustering
        new_data <- Xinput[, (clust_one[, k] == 1)]
        spear_dists <- as.matrix(stats::dist(new_data, method))
        spear_dists <- as.matrix(dist(new_data))/ (dim(new_data)[2])
        # spear_dists <- as.matrix(dist(new_data, method))
        b_vec <- c()
        #if only one element in row clust
        if(sum(indices)==1){
          a_vals <- 0
        }else{
          a_vals <- apply(spear_dists[indices, indices], 
                         1, function(x) sum(x)/(length(x)-1))
        }
        # calculate b values
        other <- (1:n_clusts_row)[-k]
        b_vals <- vector("list", length = (n_clusts_row - 1))
        t <- 1
        for(l in other){
            oth_ind <- clust_two[, l] == 1
            if(sum(oth_ind)==0){
              #if other is empty
              b_val <- rep(Inf, sum(indices))
            }else if(all(oth_ind == indices)){
              #if other is equal
              b_val <- rep(Inf, sum(indices))
            }else if((sum(oth_ind) == 1) | (sum(indices) == 1)){
              #if either has only one element
              b_val <- mean(spear_dists[indices, oth_ind])
            }else{
              b_val <- rowMeans(spear_dists[indices, oth_ind])
            }
            b_vec <- c(b_vec, mean(b_val))
            b_vals[[t]] <- b_val
            t <- t + 1
            }
          closest <- which.min(b_vec)
          b_vals <- b_vals[[closest]]
        if(all(b_vals == Inf)){
          s_vals[[k]] <- 0
          sil_score[k] <- 0
        }else if(all(b_vals == 0) & all(a_vals == 0)){
          sil_score[k] <- 0
        }else{
          s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
          s_vals[[k]] <- s_con
          sil_score[k] <- mean(s_con)
        }
    }
    }
  if (sum(sil_score) == 0){
    sil <- 0
  }else{
    sil <- ifelse(sum(sil_score != 0) == 1, sum(sil_score),
                 sum(sil_score) / (sum(sil_score != 0)) -
                 2 * sd(sil_score[sil_score != 0]))
    sil <- ifelse(sum(sil_score != 0) == 1, sum(sil_score),
                 sum(sil_score) / (sum(sil_score != 0)))
  }
  #return relationships and mean of sil_score across views
  return(list("sil" = sil, "vals" = s_vals, "repeat" = rep, "b"=b_vals, "a"=a_vals, "dists"=spear_dists, "bvec"=b_vec, "s"=s_con))
}


sil_score_inner <- function(Xinput, row_clustering, col_clustering, method="euclidean"){
    #simultaneously calculates silhouette score for a 
    #clustering as well matching clusters correctly.
    # if(any(colSums(Xinput)!=1)){
    #   Xinput <- single_alt_l1_normalisation(Xinput)$newMatrix
    # }
    n_clusts <- ncol(row_clustering)
    n_clusts_row <- n_clusts
    sil_score <- rep(0, length = n_clusts)
    clust_one <- col_clustering
    clust_two <- row_clustering
    if(n_clusts_row==1){
      clust_two <- cbind(clust_two, rbinom(nrow(row_clustering), 1, 0.1))
      n_clusts_row <- ncol(clust_two)
    }
    while(check(clust_two)){
      clust_two <- cbind(clust_two, rbinom(nrow(row_clustering), 1, 0.1))
      n_clusts_row <- ncol(clust_two)
    }
    # define whether repeat needs to happen
    rep <- ifelse(n_clusts_row == n_clusts, FALSE, TRUE)
    s_vals <- vector("list", length=n_clusts)
    for (k in 1:n_clusts){
        indices <- clust_two[, k] == 1
       # if row or col cluster empty, set to 0
        if ((sum(indices) == 0) | (sum(clust_one[,k] == 1) == 0)){
          sil_score[k] <- 0
        }else {
          #select data from specific column clustering
        new_data <- Xinput[, (clust_one[, k] == 1)]
        spear_dists <- as.matrix(stats::dist(new_data, method))
        spear_dists <- as.matrix(dist(new_data))/ (dim(new_data)[2])
        # spear_dists <- as.matrix(dist(new_data, method))
        b_vec <- c()
        #if only one element in row clust
        for()
        if(sum(indices)==1){
          a_vals <- 0
        }else{
          a_vals <- apply(spear_dists[indices, indices], 
                         1, function(x) sum(x)/(length(x)-1))
        }
        # calculate b values
        other <- (1:n_clusts_row)[-k]
        b_vals <- vector("list", length = (n_clusts_row - 1))
        t <- 1
        for(l in other){
            oth_ind <- clust_two[, l] == 1
            if(sum(oth_ind)==0){
              #if other is empty
              b_val <- rep(Inf, sum(indices))
            }else if(all(oth_ind == indices)){
              #if other is equal
              b_val <- rep(Inf, sum(indices))
            }else if((sum(oth_ind) == 1) | (sum(indices) == 1)){
              #if either has only one element
              b_val <- mean(spear_dists[indices, oth_ind])
            }else{
              b_val <- rowMeans(spear_dists[indices, oth_ind])
            }
            b_vec <- c(b_vec, mean(b_val))
            b_vals[[t]] <- b_val
            t <- t + 1
            }
          closest <- which.min(b_vec)
          b_vals <- b_vals[[closest]]
        if(all(b_vals == Inf)){
          s_vals[[k]] <- 0
          sil_score[k] <- 0
        }else if(all(b_vals == 0) & all(a_vals == 0)){
          sil_score[k] <- 0
        }else{
          s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
          s_vals[[k]] <- s_con
          sil_score[k] <- mean(s_con)
        }
    }
    }
  if (sum(sil_score) == 0){
    sil <- 0
  }else{
    sil <- ifelse(sum(sil_score != 0) == 1, sum(sil_score),
                 sum(sil_score) / (sum(sil_score != 0)) -
                 2 * sd(sil_score[sil_score != 0]))
    sil <- ifelse(sum(sil_score != 0) == 1, sum(sil_score),
                 sum(sil_score) / (sum(sil_score != 0)))
  }
  #return relationships and mean of sil_score across views
  return(list("sil" = sil, "vals" = s_vals, "repeat" = rep, "b"=b_vals, "a"=a_vals, "dists"=spear_dists, "bvec"=b_vec, "s"=s_con))
}



#old bisil score
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
  return(list("scores" = sil_score, "relations" = relations, "overall" = overall_score, "a"=a_vals, "b"=b_vals, "s"=s_mat))
}
