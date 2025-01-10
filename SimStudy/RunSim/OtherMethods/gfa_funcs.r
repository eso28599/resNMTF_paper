#apply GFA
library(GFA)
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
#function to apply GFA
gfa_apply <- function(X, k , noise = 0.5, conf = 1){
    X_in <- normalizeData(X, type = "center")[[1]]
    n_views <- length(X)
    opts <- getDefaultOpts(bicluster = TRUE)
    opts <- informativeNoisePrior(X_in, opts, noise, conf)
    gfa_res <- quiet(gfa(X_in, K = k + 5, opts))
    rows <- rep(list(apply(gfa_res$X>0, 2, as.numeric)), n_views)
    all_col <- apply(gfa_res$Z>0, 2, as.numeric)
    cols <- vector("list", length = n_views)
    #determine column results
    c1 <- 0
    for(i in 1:n_views){
        #number of columns in view i
        c2 <- c1 + dim(X_in[[i]])[2]
        cols[[i]] <- all_col[(c1 + 1):c2, ]
        c1 <- c2
    }
return(list("row_clusters" = rows, "col_clusters" = cols))
}


single_nmtf <- function(Xinput, k_min=3, k_max =8, repeats = 5, stab=TRUE){
    n_views <- length(Xinput)
    row_clusters <- vector("list", length = n_views)
    col_clusters <- vector("list", length = n_views)
    error <- c()
    for(i in 1:n_views){
      results <- restMultiNMTF_run(list(Xinput[[i]]), phi = matrix(0,1,1),
     k_min, k_max, stability=stab)
      row_clusters[[i]] <- results$row_clusters[[1]]
      col_clusters[[i]] <- results$col_clusters[[1]]
      error <- c(error, results$Error)
    }
    return(list("Error" = error,
              "row_clusters" = row_clusters,
              "col_clusters" = col_clusters))
}
