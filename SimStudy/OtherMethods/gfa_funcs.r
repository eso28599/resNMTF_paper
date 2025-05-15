# apply GFA
if (!requireNamespace("GFA", quietly = TRUE)) {
  install.packages("GFA")
}
library(GFA)
library(resnmtf)
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
# function to apply GFA
gfa_apply <- function(x, k, noise = 0.5, conf = 1) {
  x_input <- normalizeData(x, type = "center")[[1]]
  n_views <- length(x)
  opts <- getDefaultOpts(bicluster = TRUE)
  opts <- informativeNoisePrior(x_input, opts, noise, conf)
  gfa_res <- quiet(gfa(x_input, K = k + 5, opts))
  rows <- rep(list(apply(gfa_res$X > 0, 2, as.numeric)), n_views)
  all_col <- apply(gfa_res$Z > 0, 2, as.numeric)
  cols <- vector("list", length = n_views)
  # determine column results
  c1 <- 0
  for (i in 1:n_views) {
    # number of columns in view i
    c2 <- c1 + dim(x_input[[i]])[2]
    cols[[i]] <- all_col[(c1 + 1):c2, ]
    c1 <- c2
  }
  return(list("row_clusters" = rows, "col_clusters" = cols))
}


single_nmtf <- function(x_input, k_min = 3, k_max = 8,
                        repeats = 5, stab = TRUE) {
  n_views <- length(x_input)
  row_clusters <- vector("list", length = n_views)
  col_clusters <- vector("list", length = n_views)
  error <- c()
  for (i in 1:n_views) {
    results <- apply_resnmtf(list(x_input[[i]]),
      phi = matrix(0, 1, 1),
      k_min, k_max, stability = stab
    )
    row_clusters[[i]] <- results$row_clusters[[1]]
    col_clusters[[i]] <- results$col_clusters[[1]]
    error <- c(error, results$Error)
  }
  return(list(
    "Error" = error,
    "row_clusters" = row_clusters,
    "col_clusters" = col_clusters
  ))
}
