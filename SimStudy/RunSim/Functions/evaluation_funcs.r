suppressPackageStartupMessages(library(mclust))
library(clue)
library(aricode)
library(cluster)
suppressPackageStartupMessages(library(proxy))

single_alt_l1_normalisation <- function(Xmatrix){
  # normalises a matrix so that the l1 norm of each column is 1
  Q <- diag(colSums(Xmatrix))
  #solve Q
  newMatrix <- Xmatrix %*% solve(Q)
  return(list("Q" = Q, "newMatrix" = newMatrix))
}

jaccard <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <-  length(a) + length(b) - intersection
    if (union == 0 ){
      return(0)
    }else {
      return(intersection / union)
    }
}
cart_prod <- function(a, b) {
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

jaccard_res <- function(row_c, col_c, true_r, true_c, stability = FALSE){
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
      return(list("rec" = rep(0, 2), "rel" = rep(0, 2), "f_score" = rep(0, 2), "relations" = rep(0, n)))
    }
  }
  # if no biclusters present and none detected - score of 1
  if (m_0 == 0 && n_0 == 0) {
    if (stability) {
      return(1)
    }else{
      return(list("rec" = rep(1, 2), "rel" = rep(1, 2), "f_score" = rep(1, 2), "relations" = rep(0, n)))
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
  rev <- ifelse(sum(apply(jac_mat, 2, max) != 0) == 0, 0, 
    sum(apply(jac_mat, 2, max)) / n_0)
  f <- ifelse(rel * rev == 0, 0, 2 * rel * rev / (rel + rev))
  relations <- apply(jac_mat, 2, which.max)
  return(list("rec" = rep(rev, 2), "rel" = rep(rel, 2), "f_score" = rep(f, 2), "relations" = relations))
}

check <- function(matrix){

       if(sum(colSums(matrix)!=0)>1){
        matrix <- matrix[, colSums(matrix)!=0]
       }
       n_clusts <- ncol(matrix)
       equal <- diag(n_clusts)
       for(i in 1:(n_clusts-1)){
              for(j in (i+1):n_clusts){
                     check <- all(matrix[,i]==matrix[,j])
                     equal[i,j] <- check
                     equal[j,i] <- check
              }
       }
       return(nrow(unique(equal))==2)
}

# sil_score_inner <- function(Xinput, row_clustering, col_clustering, method="euclidean"){
#     #simultaneously calculates silhouette score for a 
#     #clustering as well matching clusters correctly.
#     if(any(colSums(Xinput)!=1)){
#       Xinput <- single_alt_l1_normalisation(Xinput)$newMatrix
#     }
#     n_clusts <- ncol(row_clustering)
#     n_clusts_row <- n_clusts
#     sil_score <- rep(0, length = n_clusts)
#     clust_one <- col_clustering
#     clust_two <- row_clustering
#     if(check(clust_two)){
#       clust_two <- cbind(clust_two,rbinom(nrow(row_clustering), 1, 0.1))
#       n_clusts_row <- ncol(clust_two)
#     }
#     # define whether repeat needs to happen
#     rep <- ifelse(n_clusts_row == n_clusts, FALSE, TRUE)
#     s_vals <- vector("list", length=n_clusts)
#     for (k in 1:n_clusts){
#       #select data from specific column clustering
#       new_data <- Xinput[, (clust_one[, k] == 1)]
#       spear_dists <- as.matrix(dist(new_data, method))

#       indices <- clust_two[, k] == 1
#       b_vec <- c()
#       if ((sum(indices) == 0) | (sum(clust_one[,k] == 1) == 0)){
#         sil_score[k] <- 0
#       }else {
#         if(sum(indices)==1){
#           a_vals <- 0
#         }else{
#           a_vals <- apply(spear_dists[indices, indices], 1, function(x) sum(x)/(length(x)-1))
#         }
#         #no other clusts
#         if(sum(colSums(row_clustering) != 0)==1){
#           if(sum(indices==1)){
#             b_vals <- mean(spear_dists[indices, clust_two[, k] != 1])
#           }else{
#             b_vals <- rowMeans(spear_dists[indices, clust_two[, k] != 1])
#           }
#           b_vec <- c(b_vec, mean(b_vals))
#         }else{
#           other <- (1:n_clusts_row)[-k]
#           b_vals <- vector("list", length = (n_clusts_row - 1))
#           t <- 1
#           for(l in other){
#               oth_ind <- clust_two[, l] == 1
#               if(sum(oth_ind)==0){
#                 b_val <- rep(Inf, sum(indices))
#               }else if(all(oth_ind==indices)){
#                 b_val <- rep(Inf, sum(indices))
#               }else if((sum(oth_ind)==1)|(sum(indices)==1)){
#                 b_val <- mean(spear_dists[indices, oth_ind])
#               }else{
#                 b_val <- rowMeans(spear_dists[indices, oth_ind])
#               }
#               b_vec <- c(b_vec, mean(b_val))
#               b_vals[[t]] <- b_val
#               t <- t + 1
#               }  
#           }
#           closest <- which.min(b_vec)
#           b_vals <- b_vals[[closest]]
#         if(all(b_vals==Inf)){
#           s_vals[[k]] <- 0
#           sil_score[k] <- 0
#         }else if(all(b_vals==0)&all(a_vals==0)){
#           sil_score[k] <- 0
#         }else{
#           s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
#           s_vals[[k]] <- s_con
#           sil_score[k] <- mean(s_con)
#           if(is.na(sil_score[k])){
#             print((a_vals))
#             print((b_vals))
#           }
#         }
#     }
#     }
#   if (sum(sil_score) == 0){
#     sil <- 0
#   }else{
#     sil <- ifelse(sum(sil_score != 0) == 1, sum(sil_score),
#                  sum(sil_score) / (sum(sil_score != 0)) -
#                  2 * sd(sil_score[sil_score != 0]))
#     sil <- ifelse(sum(sil_score != 0) == 1, sum(sil_score),
#                  sum(sil_score) / (sum(sil_score != 0)))
#   }
#   #return relationships and mean of sil_score across views
#   return(list("sil" = sil, "vals" = s_vals, "repeat" = rep))
# }

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
    if(check(clust_two)){
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
  # if no clusters, return 0
  if(sum(row_clustering)==0|sum(col_clustering)==0){
    return(list("sil" = 0, "vals" = 0))
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

csr <- function(row_clustering, col_clustering, true_row_clustering, true_col_clustering){
  n_row_cl <- sum(colSums(row_clustering) > 0)
  n_col_cl <- sum(colSums(col_clustering) > 0)
  true_row_cl <- sum(colSums(true_row_clustering) > 0)
  true_col_cl <- sum(colSums(true_col_clustering) > 0)

  #calculate csr
  row_csr <- 1 - abs(n_row_cl - true_row_cl) / (n_row_cl + true_row_cl + 1)
  col_csr <- 1 - abs(n_col_cl - true_col_cl) / (n_col_cl + true_col_cl + 1)
  return(c(row_csr, col_csr))
}

evaluate_simulation_comp <- function(row_clustering, col_clustering, true_row_clustering, true_col_clustering, data_views, index = 2){
  #' row_clust/col_clust: list of the clustering for each view of the rows/columns given by a method
  #' true_row/col_clustering: list of the true clustering for each view of the rows/columns 
  n_views <- length(row_clustering)
  #accuracy table for each view - 1st row is row-clustering, second is column
  csr <- matrix(0, nrow = 2, ncol = n_views)
  rownames(csr) <- c("Row-clustering", "Column-clustering")
  #relevance/recovery for each view
  rel <- matrix(0, nrow = 2, ncol = n_views)
  rownames(rel) <- c("Row-clustering", "Column-clustering")
  rec <- matrix(0, nrow = 2, ncol = n_views)
  rownames(rec) <- c("Row-clustering", "Column-clustering")
  f_score <- matrix(0, nrow = 2, ncol = n_views)
  rownames(f_score) <- c("Row-clustering", "Column-clustering")
  sil_mat <- matrix(0, nrow = 2, ncol = n_views)
  rownames(sil_mat) <- c("Row-clustering", "Column-clustering")
  relations_mat <- matrix(0, nrow = 2, ncol = n_views)
  rownames(relations_mat) <- c("Row-clustering", "Column-clustering")
  #go through each view
  #initialise storage of relations 
  relations_list <- vector("list", length = n_views)
  for (i in 1:n_views){
    sil_mat[, i] <- rep(sil_score(data_views[[i]], row_clustering[[i]],
                                  col_clustering[[i]])$sil, 2)
                    
    jaccard <- jaccard_res(row_clustering[[i]],
                                  col_clustering[[i]],
                                   true_row_clustering[[i]],
                                    true_col_clustering[[i]])
    rel[, i] <- jaccard$rel
    rec[, i] <- jaccard$rec
    f_score[, i] <- jaccard$f_score
    relations_list[[i]] <- jaccard$relations 
    #in the results
    csr[, i] <- csr(row_clustering[[i]],
                                  col_clustering[[i]],
                                   true_row_clustering[[i]],
                                    true_col_clustering[[i]])
    
  }
  #calculate relations matrix
  for(i in 1:n_views){
    score <- 0
    for(j in 1:n_views){
        score <- score + as.numeric(all(relations_list[[i]] == relations_list[[j]]))
    }
    relations_mat[, i] <- (score - 1) / max(n_views - 1, 1)
  }
  return(list("CSR" = csr,
              "BiS" = sil_mat,
              "Rel" =  rel,
              "Rec" = rec,
              "F_score" = f_score, 
              "Relations" = relations_mat))
}

eval_method <- function(data_name, file_path,  true_rows, true_cols, data_views){
    #now apply to gfa
    file_path <- paste0(data_name, file_path)
    row_filename <- paste0(file_path, "/row_clusts.xlsx")
    col_filename <- paste0(file_path, "/col_clusts.xlsx")
    results <- evaluate_simulation_comp(import_matrix(row_filename),
                                              import_matrix(col_filename),
                                              true_rows, true_cols, data_views)
    #export each data frame to separate sheets in same Excel fi;e
    path_to_save <- paste0(file_path, "_results.xlsx")
    openxlsx::write.xlsx(results, file = path_to_save)
}

jaccard_row<- function(row_c, true_r, print=FALSE){
  m <- ncol(row_c)
  n <- ncol(true_r)
  # if no biclusters detected but some are present
  # return 0
  # if no biclusters present but some are detected - score of 0
  m_0 <- sum(colSums(row_c) != 0) #no of clusters actually detected
  n_0 <- sum(colSums(true_r) != 0) #no of true clusters
  if ((m_0 == 0 && n_0 != 0) || (n_0 == 0 && m_0 != 0)) {
    return(list("rec" = rep(0, 2), "rel" = rep(0, 2), "f_score" = rep(0, 2), "relations" = rep(0, n)))
  }
  # if no biclusters present and none detected - score of 1
  if (m_0 == 0 && n_0 == 0) {
    return(list("rec" = rep(1, 2), "rel" = rep(1, 2), "f_score" = rep(1, 2), "relations" = rep(0, n)))
  }
  rows <- 1:nrow(row_c)
  #initialise storage of jaccard index between pairs
  jac_mat <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m){
    r_i <- rows[row_c[, i] == 1]
    for (j in 1:n){
        tr_i <- rows[true_r[, j] == 1]
        jac_mat[i, j] <- jaccard(r_i, tr_i)
    }
  }
  if(print){
    print(jac_mat)
    print(apply(jac_mat, 1, max))
    print(apply(jac_mat, 2, max))
  }
  rel <- ifelse(sum(apply(jac_mat, 1, max) != 0) == 0, 0,
    sum(apply(jac_mat, 1, max)) / m_0)
  rev <- ifelse(sum(apply(jac_mat, 2, max) != 0) == 0, 0, 
    sum(apply(jac_mat, 2, max)) / n_0)
  f <- ifelse(rel * rev == 0, 0, 2 * rel * rev / (rel + rev))
  relations <- apply(jac_mat, 2, which.max)
  return(list("rec" = rep(rev, 2), "rel" = rep(rel, 2), "f_score" = rep(f, 2), "relations" = relations))
}

get_sil_mean <- function(vec){
  return(ifelse(sum(vec)==0, 0, mean(vec[vec!=0])))
}

calc_all_sils <- function(data,res){
      n_views <- length(data)
      bisils_euc <- c()
      bisils_man <- c()
      bisils_cosine <- c()
      for(i in 1:n_views){
        bisils_euc <- c(bisils_euc, 
                        sil_score(data[[i]],res$row_clusters[[i]],res$col_clusters[[i]])$sil)
        bisils_man <- c(bisils_man, 
                        sil_score(data[[i]],res$row_clusters[[i]],res$col_clusters[[i]], method="manhattan")$sil)
        bisils_cosine <- c(bisils_cosine,
                    sil_score(data[[i]],res$row_clusters[[i]],res$col_clusters[[i]], method="cosine")$sil)
      }
      bisils_euc <- c(bisils_euc, get_sil_mean(bisils_euc))
      bisils_man <- c(bisils_man, get_sil_mean(bisils_man))
      bisils_cosine <- c(bisils_cosine, get_sil_mean(bisils_cosine))
      return(list("euc"=bisils_euc, "man"=bisils_man, "cosine"=bisils_cosine))
}

all_jaccs <- function(rows, res){
  n_views <- length(res)
  jaccs <- vector("list", length=n_views)
  for(i in 1:n_views){
    jaccs[[i]] <- jaccard_row(res[[i]],rows)
  }
  #scores
  fs <- c()
  rels <- c()
  recs <- c()
  for(i in 1:n_views){
    fs <- c(fs, jaccs[[i]]$f_score[1])
    rels <- c(rels, jaccs[[i]]$rel[1])
    recs <- c(recs, jaccs[[i]]$rec[1])
  }
  return(c(fs, mean(fs), rels, mean(rels), recs, mean(recs)))
}


dis_results <- function(data, rows, res, phi ,rep, path, row_same=FALSE){
        #save data
        openxlsx::write.xlsx(res$row_clusters,
                file = paste0(path, "/data/row_clusts", phi, "_", rep, ".xlsx"))
        openxlsx::write.xlsx(res$col_clusters,
                file = paste0(path, "/data/col_clusts", phi, "_", rep, ".xlsx"))
        sils <- calc_all_sils(data, res)
        #jaccard
        if(row_same){
          jaccs <- all_jaccs(rows, res$row_clusters)
        }else{
          jaccs <- all_jaccs(rows, res$col_clusters)
        }
        
        #no clusts
        no_clusts <- max(sapply(res$row_clusters,
                                 function(x) sum(colSums(x) != 0)))
        return(c(jaccs, sils$euc, sils$cosine, sils$man, no_clusts))
}