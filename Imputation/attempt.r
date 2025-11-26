source("SimStudy/Functions/data_generation.r")
library(rje)
library(hash)
n_views <- 2
data <- multi_view(c(10, 10), c(10, 10), 3, 2, 5, seed = 1)
# list of the row names for each view
row_names <- list(paste0("row_", 1:10), paste0("row_", 5:14))
col_names <- list(paste0("col_", 1:10), paste0("col_", 5:14))
for (view in 1:n_views) {
  rownames(data$data_views[[view]]) <- row_names[[view]]
  colnames(data$data_views[[view]]) <- col_names[[view]]
}

produce_indices <- function(
    row_lists, col_lists,
    existing_view_subsets_r, existing_view_subsets_c, n_views) {
  shared_rows <- vector("list", length = n_views)
  shared_cols <- vector("list", length = n_views)
  for (view1 in 1:n_views) {
    shared_v1 <- hash()
    shared_c_v1 <- hash()
    for (view2 in (1:3)[-view1]) {
      common_rows <- unlist(
        row_lists[
          sapply(
            existing_view_subsets_r,
            function(x) all(c(view1, view2) %in% x)
          )
        ]
      )
      common_cols <- unlist(
        col_lists[
          sapply(
            existing_view_subsets_c,
            function(x) all(c(view1, view2) %in% x)
          )
        ]
      )
      if (length(common_rows) != 0) {
        shared_v1[[as.character(view2)]] <- common_rows
      } else {
        shared_v1[[as.character(view2)]] <- NA
      }
      if (length(common_cols) != 0) {
        shared_c_v1[[as.character(view2)]] <- common_rows
      } else {
        shared_c_v1[[as.character(view2)]] <- NA
      }
    }
    shared_rows[[view1]] <- shared_v1
    shared_cols[[view1]] <- shared_c_v1
  }
  return(list("shared_rows" = shared_rows, "sahred_cols" = shared_cols))
}

reorder_data <- function(data, n_views) {
  # for each view v, find R^(v) - the list of existing view subsets A for view v
  # define relevant power set, not including the empty subset
  # a row is in existing view subset A
  row_names <- lapply(data, rownames)
  col_names <- lapply(data, colnames)
  power_set <- powerSetCond(1:n_views)
  existing_view_subsets_r <- list()
  existing_view_subsets_c <- list()
  row_lists <- list()
  col_lists <- list()
  # extract relevant row_names
  for (views in power_set) {
    # define {1, ..., n} \ A
    negate_v_sub <- (1:n_views)[!((1:n_views) %in% views)]
    # find rows belonging to views v in A
    rows <- Reduce(intersect, row_names[views])
    # find rows belonging to views v not in A
    rows_negation <- Reduce(union, row_names[negate_v_sub])
    # select the rows present in A but not in {1, ..., n} \ A
    rows_in_a <- setdiff(rows, rows_negation)
    print(negate_v_sub)
    print(rows)
    print(rows_negation)
    if (length(rows_in_a) != 0) {
      existing_view_subsets_r <- c(existing_view_subsets_r, list(views))
      row_lists <- c(row_lists, list(rows_in_a))
    }
    # find cols belonging to views v in A
    cols <- Reduce(intersect, col_names[views])
    # find cols belonging to views v not in A
    cols_negation <- Reduce(union, col_names[negate_v_sub])
    # select the cols present in A but not in {1, ..., n} \ A
    cols_in_a <- setdiff(cols, cols_negation)
    if (length(cols_in_a) != 0) {
      existing_view_subsets_c <- c(existing_view_subsets_c, list(views))
      col_lists <- c(col_lists, list(cols_in_a))
    }
  }

  # reorder_data
  data_reordered <- vector("list", length = n_views)
  row_orders <- list()
  col_orders <- list()
  for (view in 1:n_views) {
    row_order <- unlist(
      row_lists[sapply(existing_view_subsets_r, function(x) view %in% x)]
    )
    col_order <- unlist(
      col_lists[sapply(existing_view_subsets_c, function(x) view %in% x)]
    )
    data_reordered[[view]] <- data[[view]][row_order, col_order]
    row_orders <- c(row_orders, list(row_order))
    col_orders <- c(col_orders, list(col_order))
  }

  indices <- produce_indices(
    row_lists, col_lists,
    existing_view_subsets_r, existing_view_subsets_c, n_views
  )
  return(list(
    "data_reordered" = data_reordered,
    "existing_view_subsets_r" = existing_view_subsets_r,
    "existing_view_subsets_c" = existing_view_subsets_c,
    "row_indices" = indices$shared_rows,
    "col_indices" = indices$shared_cols,
    "row_lists" = row_lists,
    "col_lists" = col_lists
  ))
}

data <- multi_view(c(10, 10), c(10, 10), 3, 2, 5, seed = 1)
# list of the row names for each view
row_names <- list(paste0("row_", 1:10), paste0("row_", 1:10))
col_names <- list(paste0("col_", 1:10), paste0("col_", 5:14))
for (view in 1:n_views) {
  rownames(data$data_views[[view]]) <- row_names[[view]]
  colnames(data$data_views[[view]]) <- col_names[[view]]
}
reordered_list <- reorder_data(data$data_views, n_views)

# tests for data reordering
if (
  length(Reduce(intersect, row_lists)) != 0 |
    (!setequal(Reduce(union, row_lists), Reduce(union, row_names)))
) {
  print("Row partition not implemented correctly.")
}
if (
  length(Reduce(intersect, col_lists)) != 0 |
    (!setequal(Reduce(union, col_lists), Reduce(union, col_names)))
) {
  print("Column partition not implemented correctly.")
}
if (any(sapply(row_orders, length) != sapply(data$data_views, nrow))) {
  print("Row reordering not implemented correctly.")
}
if (any(sapply(col_orders, length) != sapply(data$data_views, ncol))) {
  print("Column reordering not implemented correctly.")
}
for (view in 1:n_views) {
  if (!setequal(
    rownames(data$data_views[[views]]),
    rownames(data_reordered[[views]])
  )) {
    print("Row reordering not implemented correctly.")
  }
  if (!setequal(
    colnames(data$data_views[[views]]),
    colnames(data_reordered[[views]])
  )) {
    print("Column reordering not implemented correctly.")
  }
}

give_names <- function(data) {
  names_missing_rows <- sapply(data, function(x) is.null(rownames(x)))
  names_missing_cols <- sapply(data, function(x) is.null(colnames(x)))
  # if all rows aren't named - name
  if (all(names_missing_rows)) {
    n <- 1
    for (i in 1:n_views) {
      if (is.null(rownames(data[[i]]))) {
        rownames(data[[i]]) <- paste0("row_", n:(n - 1 + nrow(data[[i]])))
        n <- n + nrow(data[[i]])
      }
    }
    # if some views are completely missing names, user must specify
  } else if (any(names_missing_rows)) {
    stop("At least one view is missing row names. Please name missing rows.")
  } else {
    # if some rows are missing names within a view but others are provided
    # user must specify
    if (any(sapply(data, function(x) any(is.na(rownames(x)))))) {
      stop("Some rows missing names. Check row names.")
    }
  }
  # now check columns
  if (all(names_missing_cols)) {
    p <- 1
    for (i in 1:n_views) {
      if (is.null(colnames(data[[i]]))) {
        colnames(data[[i]]) <- paste0("col_", p:(p - 1 + ncol(data[[i]])))
        p <- p + ncol(data[[i]])
      }
    }
  } else if (any(names_missing_cols)) {
    stop(
      "At least one view is missing column names. Please name missing columns."
    )
  } else {
    # if any columns are missing names
    if (any(sapply(data, function(x) any(is.na(colnames(x)))))) {
      stop("Some columns missing names. Check columns names.")
    }
  }
  return(list(
    "data" = data,
    "row_names" = lapply(data, rownames),
    "col_names" = lapply(data, colnames)
  ))
}

# -----------------------------
# tests for naming functions
# -----------------------------

# one row missing a name
data <- multi_view(c(10, 10), c(10, 10), 3, 2, 5, seed = 1)$data_views
data <- list(
  abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10)))),
  abs(MASS::mvrnorm(10, mu = rep(0, 10), Sigma = diag(rep(1, 10))))
)

og_row_names <- lapply(data, rownames)
og_col_names <- lapply(data, colnames)
data <- give_names(data)$data
reordering <- reorder_data(data, 2)


rownames(data[[1]]) <- c(paste0("row_", 1:9), NA)
rownames(data[[2]]) <- paste0("row_", 5:14)
colnames(data[[1]]) <- paste0("col_", 1:10)
colnames(data[[2]]) <- paste0("col_", 5:14)
expect_error(
  give_names(data),
  "Some rows missing names. Check row names."
)

# one view rows not named
data <- multi_view(c(10, 10), c(10, 10), 3, 2, 5, seed = 1)$data_views
rownames(data[[1]]) <- paste0("row_", 1:10)
colnames(data[[1]]) <- paste0("col_", 1:10)
colnames(data[[2]]) <- paste0("col_", 5:14)
expect_error(
  give_names(data),
  "At least one view is missing row names. Please name missing rows."
)

# one column missing a name
# "Some columns missing names. Check columns names."
data <- multi_view(c(10, 10), c(10, 10), 3, 2, 5, seed = 1)$data_views
rownames(data[[1]]) <- paste0("row_", 1:10)
rownames(data[[2]]) <- paste0("row_", 5:14)
colnames(data[[1]]) <- c(paste0("col_", 1:9), NA)
colnames(data[[2]]) <- paste0("col_", 5:14)
expect_error(
  give_names(data),
  "Some columns missing names. Check columns names."
)

# one view columns not named
data <- multi_view(c(10, 10), c(10, 10), 3, 2, 5, seed = 1)$data_views
rownames(data[[1]]) <- paste0("row_", 1:10)
rownames(data[[2]]) <- paste0("row_", 5:14)
colnames(data[[1]]) <- paste0("col_", 1:10)
expect_error(
  give_names(data),
  "At least one view is missing column names. Please name missing columns."
)

# no row names
data <- multi_view(c(10, 10), c(10, 10), 3, 2, 5, seed = 1)$data_views
colnames(data[[1]]) <- paste0("col_", 1:10)
colnames(data[[2]]) <- paste0("col_", 5:14)
res <- give_names(data)
expect_equal(res$row_names[[1]], paste0("row_", 1:10))
expect_equal(res$row_names[[2]], paste0("row_", 11:20))

# no column names
data <- multi_view(c(10, 10), c(10, 10), 3, 2, 5, seed = 1)$data_views
rownames(data[[1]]) <- paste0("row_", 1:10)
rownames(data[[2]]) <- paste0("row_", 5:14)
res <- give_names(data)
expect_equal(res$col_names[[1]], paste0("col_", 1:10))
expect_equal(res$col_names[[2]], paste0("col_", 11:20))

# all named
data <- multi_view(c(10, 10), c(10, 10), 3, 2, 5, seed = 1)$data_views
rownames(data[[1]]) <- paste0("row_", 1:10)
rownames(data[[2]]) <- paste0("row_", 5:14)
colnames(data[[1]]) <- paste0("col_", 1:10)
colnames(data[[2]]) <- paste0("col_", 5:14)
res <- give_names(data)
expect_equal(res$row_names[[1]], paste0("row_", 1:10))
expect_equal(res$row_names[[2]], paste0("row_", 5:14))
expect_equal(res$col_names[[1]], paste0("col_", 1:10))
expect_equal(res$col_names[[2]], paste0("col_", 5:14))
