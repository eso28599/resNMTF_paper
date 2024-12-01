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
        spear_dists <- as.matrix(dist(new_data, method))
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