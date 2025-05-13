suppressPackageStartupMessages(library(mclust))
library(clue)
library(aricode)
library(cluster)
suppressPackageStartupMessages(library(proxy))


jaccard <- function(a, b) {
  # Arguments:
  #  a: vector of elements
  #  b: vector of elements
  # Returns:
  #  Jaccard index between a and b
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  if (union == 0) {
    return(0)
  } else {
    return(intersection / union)
  }
}
cart_prod <- function(a, b) {
  # Arguments:
  #  a: vector of elements
  #  b: vector of elements
  # Returns:
  #  Cartesian product of a and b as a vector of strings
  prod <- c()
  # check a or b are not empty sets
  if (length(a) == 0 || length(b) == 0) {
    return(NULL)
  } else {
    for (k in seq_len(length(a))) {
      prod <- c(prod, paste(a[k], b))
    }
    return(prod)
  }
}

jaccard_res <- function(row_c, col_c, true_r, true_c) {
  # Arguments:
  #  row_c: clustering of rows
  #  col_c: clustering of columns
  #  true_r: true clustering of rows
  #  true_c: true clustering of columns
  # Returns:
  #  F score, recovery and relevance between found and true biclusters
  m <- ncol(row_c)
  n <- ncol(true_r)
  m_0 <- sum(colSums(row_c) != 0) # no of clusters actually detected
  n_0 <- sum(colSums(true_r) != 0) # no of true clusters

  # if no biclusters detected but some are present or
  # if no biclusters present but some are detected
  # return 0
  if ((m_0 == 0 && n_0 != 0) || (n_0 == 0 && m_0 != 0)) {
    return(list(
      "rec" = rep(0, 2),
      "rel" = rep(0, 2),
      "f_score" = rep(0, 2),
      "relations" = rep(0, n)
    ))
  }
  # if no biclusters present and none detected - score of 1
  if (m_0 == 0 && n_0 == 0) {
    return(list(
      "rec" = rep(1, 2),
      "rel" = rep(1, 2),
      "f_score" = rep(1, 2),
      "relations" = rep(0, n)
    ))
  }
  samps <- seq_len(nrow(row_c))
  feats <- seq_len(nrow(col_c))
  # initialise storage of jaccard index between pairs
  jac_mat <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    r_i <- samps[row_c[, i] == 1]
    c_i <- feats[col_c[, i] == 1]
    m_i <- cart_prod(r_i, c_i)
    for (j in 1:n) {
      tr_i <- samps[true_r[, j] == 1]
      tc_i <- feats[true_c[, j] == 1]
      m_j <- cart_prod(tr_i, tc_i)
      jac_mat[i, j] <- jaccard(m_i, m_j)
    }
  }
  rel <- ifelse(sum(apply(jac_mat, 1, max) != 0) == 0, 0,
    sum(apply(jac_mat, 1, max)) / m_0
  )
  rev <- ifelse(sum(apply(jac_mat, 2, max) != 0) == 0, 0,
    sum(apply(jac_mat, 2, max)) / n_0
  )
  f <- ifelse(rel * rev == 0, 0, 2 * rel * rev / (rel + rev))
  relations <- apply(jac_mat, 2, which.max)
  return(list(
    "rec" = rep(rev, 2),
    "rel" = rep(rel, 2),
    "f_score" = rep(f, 2),
    "relations" = relations
  ))
}


csr <- function(row_clustering, col_clustering,
                true_row_clustering, true_col_clustering) {
  n_row_cl <- sum(colSums(row_clustering) > 0)
  n_col_cl <- sum(colSums(col_clustering) > 0)
  true_row_cl <- sum(colSums(true_row_clustering) > 0)
  true_col_cl <- sum(colSums(true_col_clustering) > 0)

  # calculate csr
  row_csr <- 1 - abs(n_row_cl - true_row_cl) / (n_row_cl + true_row_cl + 1)
  col_csr <- 1 - abs(n_col_cl - true_col_cl) / (n_col_cl + true_col_cl + 1)
  return(c(row_csr, col_csr))
}

evaluate_simulation_comp <- function(row_clustering, col_clustering,
                                     true_row_clustering,
                                     true_col_clustering,
                                     data_views) {
  # Arguments:
  #  row_clustering: clustering of rows
  #  col_clustering: clustering of columns
  #  true_row_clustering: true clustering of rows
  #  true_col_clustering: true clustering of columns
  #  data_views: list of data views
  # Returns:
  #  CSR, relevance, recovery and F score between found and true biclusters
  #  for each view
  #  bisilhouette for each view
  #  relations matrix for each view
  n_views <- length(row_clustering)
  # accuracy table for each view - 1st row is row-clustering, second is column
  csr <- matrix(0, nrow = 2, ncol = n_views)
  rownames(csr) <- c("Row-clustering", "Column-clustering")
  # relevance/recovery for each view
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
  # go through each view
  # initialise storage of relations
  relations_list <- vector("list", length = n_views)
  for (i in 1:n_views) {
    sil_mat[, i] <- rep(bisilhouette(
      data_views[[i]], row_clustering[[i]],
      col_clustering[[i]]
    )$bisil, 2)

    jaccard <- jaccard_res(
      row_clustering[[i]],
      col_clustering[[i]],
      true_row_clustering[[i]],
      true_col_clustering[[i]]
    )
    rel[, i] <- jaccard$rel
    rec[, i] <- jaccard$rec
    f_score[, i] <- jaccard$f_score
    relations_list[[i]] <- jaccard$relations
    # in the results
    csr[, i] <- csr(
      row_clustering[[i]],
      col_clustering[[i]],
      true_row_clustering[[i]],
      true_col_clustering[[i]]
    )
  }
  # calculate relations matrix
  for (i in 1:n_views) {
    score <- 0
    for (j in 1:n_views) {
      score <- score +
        as.numeric(all(relations_list[[i]] == relations_list[[j]]))
    }
    relations_mat[, i] <- (score - 1) / max(n_views - 1, 1)
  }
  return(list(
    "CSR" = csr,
    "BiS" = sil_mat,
    "Rel" = rel,
    "Rec" = rec,
    "F_score" = f_score,
    "Relations" = relations_mat
  ))
}

eval_method <- function(data_name, file_path,
                        true_rows, true_cols,
                        data_views) {
  # now apply to gfa
  file_path <- paste0(data_name, file_path)
  row_filename <- paste0(file_path, "/row_clusts.xlsx")
  col_filename <- paste0(file_path, "/col_clusts.xlsx")
  results <- evaluate_simulation_comp(
    import_matrix(row_filename),
    import_matrix(col_filename),
    true_rows, true_cols, data_views
  )
  # export each data frame to separate sheets in same Excel fi;e
  path_to_save <- paste0(file_path, "_results.xlsx")
  openxlsx::write.xlsx(results, file = path_to_save)
}

jaccard_row <- function(row_c, true_r, print = FALSE) {
  m <- ncol(row_c)
  n <- ncol(true_r)
  # if no biclusters detected but some are present
  # return 0
  # if no biclusters present but some are detected - score of 0
  m_0 <- sum(colSums(row_c) != 0) # no of clusters actually detected
  n_0 <- sum(colSums(true_r) != 0) # no of true clusters
  if ((m_0 == 0 && n_0 != 0) || (n_0 == 0 && m_0 != 0)) {
    return(list(
      "rec" = rep(0, 2),
      "rel" = rep(0, 2),
      "f_score" = rep(0, 2),
      "relations" = rep(0, n)
    ))
  }
  # if no biclusters present and none detected - score of 1
  if (m_0 == 0 && n_0 == 0) {
    return(list(
      "rec" = rep(1, 2),
      "rel" = rep(1, 2),
      "f_score" = rep(1, 2),
      "relations" = rep(0, n)
    ))
  }
  rows <- seq_len(row_c)
  # initialise storage of jaccard index between pairs
  jac_mat <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m) {
    r_i <- rows[row_c[, i] == 1]
    for (j in 1:n) {
      tr_i <- rows[true_r[, j] == 1]
      jac_mat[i, j] <- jaccard(r_i, tr_i)
    }
  }
  if (print) {
    print(jac_mat)
    print(apply(jac_mat, 1, max))
    print(apply(jac_mat, 2, max))
  }
  rel <- ifelse(sum(apply(jac_mat, 1, max) != 0) == 0, 0,
    sum(apply(jac_mat, 1, max)) / m_0
  )
  rev <- ifelse(sum(apply(jac_mat, 2, max) != 0) == 0, 0,
    sum(apply(jac_mat, 2, max)) / n_0
  )
  f <- ifelse(rel * rev == 0, 0, 2 * rel * rev / (rel + rev))
  relations <- apply(jac_mat, 2, which.max)
  return(list(
    "rec" = rep(rev, 2),
    "rel" = rep(rel, 2),
    "f_score" = rep(f, 2),
    "relations" = relations
  ))
}

get_sil_mean <- function(vec) {
  return(ifelse(sum(vec) == 0, 0, mean(vec[vec != 0])))
}

calc_all_sils <- function(data, res) {
  n_views <- length(data)
  bisils_euc <- c()
  bisils_man <- c()
  bisils_cosine <- c()
  for (i in 1:n_views) {
    bisils_euc <- c(
      bisils_euc,
      bisilhouette(
        data[[i]],
        res$row_clusters[[i]],
        res$col_clusters[[i]]
      )$bisil
    )
    bisils_man <- c(
      bisils_man,
      bisilhouette(data[[i]],
        res$row_clusters[[i]],
        res$col_clusters[[i]],
        method = "manhattan"
      )$bisil
    )
    bisils_cosine <- c(
      bisils_cosine,
      bisilhouette(data[[i]],
        res$row_clusters[[i]],
        res$col_clusters[[i]],
        method = "cosine"
      )$bisil
    )
  }
  bisils_euc <- c(bisils_euc, get_sil_mean(bisils_euc))
  bisils_man <- c(bisils_man, get_sil_mean(bisils_man))
  bisils_cosine <- c(bisils_cosine, get_sil_mean(bisils_cosine))
  return(list("euc" = bisils_euc, "man" = bisils_man, "cosine" = bisils_cosine))
}

all_jaccs <- function(rows, res) {
  n_views <- length(res)
  jaccs <- vector("list", length = n_views)
  for (i in 1:n_views) {
    jaccs[[i]] <- jaccard_row(res[[i]], rows)
  }
  # scores
  fs <- c()
  rels <- c()
  recs <- c()
  for (i in 1:n_views) {
    fs <- c(fs, jaccs[[i]]$f_score[1])
    rels <- c(rels, jaccs[[i]]$rel[1])
    recs <- c(recs, jaccs[[i]]$rec[1])
  }
  return(c(fs, mean(fs), rels, mean(rels), recs, mean(recs)))
}


dis_results <- function(data, rows, res, phi, rep, path, row_same = FALSE) {
  # save data
  openxlsx::write.xlsx(res$row_clusters,
    file = paste0(path, "/data/row_clusts", phi, "_", rep, ".xlsx")
  )
  openxlsx::write.xlsx(res$col_clusters,
    file = paste0(path, "/data/col_clusts", phi, "_", rep, ".xlsx")
  )
  sils <- calc_all_sils(data, res)
  # jaccard
  if (row_same) {
    jaccs <- all_jaccs(rows, res$row_clusters)
  } else {
    jaccs <- all_jaccs(rows, res$col_clusters)
  }

  # no clusts
  no_clusts <- max(sapply(
    res$row_clusters,
    function(x) sum(colSums(x) != 0)
  ))
  return(c(jaccs, sils$euc, sils$cosine, sils$man, no_clusts))
}
