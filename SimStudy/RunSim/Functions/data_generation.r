library(openxlsx)
library(MASS)
library(Matrix)

# toy examples
get_dims <- function(n_r, n_c, k, row_e = 1, col_e = 1, row_o = 0, col_o = 0) {
    #' row_e portion of rows in biclusters
    part_r <- floor((row_e * n_r) / 7)
    part_c <- floor((col_e * n_c) / 7)
    if (row_e == 1) {
        if (k == 2) {
            rows <- c(4 * part_r, n_r - 4 * part_r)
        } else if (k == 3) {
            rows <- c(3 * part_r, 2 * part_r, n_r - 5 * part_r)
        } else if (k == 4) {
            rows <- c(3 * part_r, 2 * part_r, part_r, n_r - 6 * part_r)
        } else if (k == 5) {
            rows <- c(2 * part_r, 2 * part_r, part_r, part_r, n_r - 6 * part_r)
        } else if (k == 6) {
            rows <- c(2 * part_r, part_r, part_r, part_r, part_r, n_r - 6 * part_r)
        }
    } else {
        if (k == 2) {
            rows <- c(4 * part_r, 3 * part_r)
        } else if (k == 3) {
            rows <- c(3 * part_r, 2 * part_r, 2 * part_r)
        } else if (k == 4) {
            rows <- c(3 * part_r, 2 * part_r, part_r, part_r)
        } else if (k == 5) {
            rows <- c(2 * part_r, 2 * part_r, part_r, part_r, part_r)
        } else if (k == 6) {
            rows <- c(2 * part_r, part_r, part_r, part_r, part_r, part_r)
        }
    }

    if (col_e == 1) {
        if (k == 2) {
            cols <- c(4 * part_c, n_c - 4 * part_c)
        } else if (k == 3) {
            cols <- c(3 * part_c, 2 * part_c, n_c - 5 * part_c)
        } else if (k == 4) {
            cols <- c(3 * part_c, 2 * part_c, part_c, n_c - 6 * part_c)
        } else if (k == 5) {
            cols <- c(2 * part_c, 2 * part_c, part_c, part_c, n_c - 6 * part_c)
        } else if (k == 6) {
            cols <- c(2 * part_c, part_c, part_c, part_c, part_c, n_c - 6 * part_c)
        }
    } else {
        if (k == 2) {
            cols <- c(4 * part_c, 3 * part_c)
        } else if (k == 3) {
            cols <- c(3 * part_c, 2 * part_c, 2 * part_c)
        } else if (k == 4) {
            cols <- c(3 * part_c, 2 * part_c, part_c, part_c)
        } else if (k == 5) {
            cols <- c(2 * part_c, 2 * part_c, part_c, part_c, part_c)
        } else if (k == 6) {
            cols <- c(2 * part_c, part_c, part_c, part_c, part_c, part_c)
        }
    }


    rows <- sort(rows, decreasing = TRUE)
    cols <- sort(cols, decreasing = TRUE)

    row_overlap <- floor(rows * row_o)
    row_overlap[k] <- 0
    col_overlap <- floor(cols * col_o)
    col_overlap[k] <- 0

    row_start <- c(1, cumsum(rows)[-k] + 1)
    if (part_r == 0) {
        row_start <- rep(0, k)
    }
    row_end <- cumsum(rows)
    row_end <- row_end + row_overlap
    col_start <- c(1, cumsum(cols)[-k] + 1)
    if (part_c == 0) {
        col_start <- rep(0, k)
    }
    col_end <- cumsum(cols)
    col_end <- col_end + col_overlap

    return(list(
        "Row_s" = row_start, "Row_e" = row_end,
        "Col_s" = col_start, "Col_e" = col_end
    ))
}

one_view_adv <- function(row_dims, col_dims, k, noise, signal, row_e = 1, col_e = 1, row_o = 0, col_o = 0) {
    #' row_dims: vector of sizes of each row cluster in this view
    #' col_dims: vector of sizes of each row cluster in this view
    #' noise: variance of the noise added to views
    #'
    #' view: matrix of this dataview
    #' truth_row: vector indicating membership of row clusters
    #' truth_col: vector indicating membership of col clusters

    # generate noise
    X_noise <- mvrnorm(n = row_dims, mu = rep(0, col_dims), Sigma = diag(noise, col_dims))

    # Create a list as input for block-diagonal data generation
    # list length of no of clusters in this first view

    # generate a mvn with n=number of individuals of view and no of features
    # equal to number of features of the view
    # mean of 5 for each column
    # covariance matrix is identity - each feature is independent
    # define vectors indicating true row/column cluster membership for each view
    true_row <- matrix(0, nrow = row_dims, ncol = k)
    true_col <- matrix(0, nrow = col_dims, ncol = k)
    X_view <- matrix(0, nrow = row_dims, ncol = col_dims)
    dims <- get_dims(row_dims, col_dims, k, row_e, col_e, row_o, col_o)
    row_start <- dims$Row_s
    row_end <- dims$Row_e
    col_start <- dims$Col_s
    col_end <- dims$Col_e
    # row_start
    for (i in 1:k) {
        n_r <- (row_end[i] - row_start[i] + 1)
        n_c <- (col_end[i] - col_start[i] + 1)
        if (n_r != 0) {
            X_view[(row_start[i]):(row_end[i]), (col_start[i]):(col_end[i])] <- mvrnorm(n = n_r, mu = rep(signal, n_c), Sigma = diag(n_c))
            # X_view[(row_start[i]):(row_end[i]), (col_start[i]):(col_end[i])] <-  mvrnorm(n = n_r, mu = rep(100, n_c), Sigma = diag(n_c))

            true_row[(row_start[i]):(row_end[i]), i] <- 1
            true_col[(col_start[i]):(col_end[i]), i] <- 1
        }
    }
    # add noise to first view
    X_view <- abs(X_view) + abs(X_noise)
    return(list(view = X_view, truth_row = true_row, truth_col = true_col))
}

make_longer <- function(vec, n) {
    if (length(vec) == 1) {
        vec <- rep(vec, n)
    }
    return(vec)
}

multi_view <- function(row_dims, col_dims, k, noise, signal, row_e = 1, col_e = 1, row_o = 0, col_o = 0, row_same_shuffle = TRUE, col_same_shuffle = TRUE, seed = FALSE, file_path = NA) {
    #' rowClusters: n length list of vectors of row cluster sizes in each view
    #' rowClusters: n length list of vectors of column cluster sizes in each view
    #' seed: logical indicator, default is FALSE, if true sets seed so same data is generated each time
    #'
    #' data_views: n length list of data views
    #' truth_rows: n length list of vectors of true row cluster membership for each view
    #' col_rows: n length list of vectors of true column cluster membership for each view

    if (is.numeric(seed)) {
        set.seed(seed)
    }
    n_views <- length(row_dims)
    # Introduce simulated data-views- where we store the views
    X_trial <- vector("list", length = n_views)
    # list to store row clusterings for each dataview
    true_row_clusterings <- vector("list", length = n_views)
    true_col_clusterings <- vector("list", length = n_views)
    # initialise index to shuffle rows and columns
    # new row/col index
    # if shuffle is true - dims of each view must be the same and represent same objects
    if (row_same_shuffle) {
        new_row_ind <- sample(row_dims[1])
    }
    if (col_same_shuffle) {
        new_col_ind <- sample(col_dims[1])
    }
    row_e <- make_longer(row_e, n_views)
    row_o <- make_longer(row_o, n_views)
    col_e <- make_longer(col_e, n_views)
    col_o <- make_longer(col_o, n_views)
    # generate data and store
    for (i in 1:n_views) {
        if (!row_same_shuffle) {
            new_row_ind <- sample(row_dims[i])
        }
        if (!col_same_shuffle) {
            new_col_ind <- sample(col_dims[i])
        }
        data_i <- one_view_adv(
            row_dims[i], col_dims[i],
            k, noise, signal, row_e[i], col_e[i], row_o[i], col_o[i]
        )
        if (!is.na(file_path)) {
            pdf(paste0(file_path, "/true_image_", i, ".pdf"))
            image(t(data_i$view))
            dev.off()
        }
        X_trial[[i]] <- (data_i$view)[new_row_ind, new_col_ind]
        true_row_clusterings[[i]] <- (data_i$truth_row)[new_row_ind, ]
        true_col_clusterings[[i]] <- (data_i$truth_col)[new_col_ind, ]
    }

    return(list(data_views = X_trial, truth_rows = true_row_clusterings, truth_cols = true_col_clusterings))
}

save_data <- function(row_dims, col_dims, k, file_path, noise, row_e = 1, col_e = 1, row_o = 0, col_o = 0, row_same_shuffle = TRUE, col_same_shuffle = TRUE, signal = 5) {
    # can change the noise parameter here for level of noise in views
    data <- multi_view(row_dims, col_dims, k, noise, signal, row_e, col_e, row_o, col_o, row_same_shuffle, col_same_shuffle)
    # save data as a file in given directory
    # export each data frame to separate sheets in same Excel file
    openxlsx::write.xlsx(data$data_views, file = paste0(file_path, "/data.xlsx")) # nolint
    openxlsx::write.xlsx(data$truth_rows, file = paste0(file_path, "/true_rows.xlsx")) # nolint: line_length_linter.
    openxlsx::write.xlsx(data$truth_cols, file = paste0(file_path, "/true_cols.xlsx")) # nolint: line_length_linter.
}

save_data_noise <- function(row_dims, col_dims, k, file_path, method_vec, noise_vec, row_same_shuffle = TRUE, col_same_shuffle = TRUE) {
    # can change the noise parameter here for level of noise in views
    data <- multi_view(row_dims, col_dims, k, 0, 5, row_same_shuffle = row_same_shuffle, col_same_shuffle = col_same_shuffle)
    # save data as a file in given directory
    # export each data frame to separate sheets in same Excel file
    for (i in 1:length(noise_vec)) {
        file_path2 <- paste0(file_path, paste0("/", method_vec[i]))
        data_set <- vector("list", length = length(row_dims))
        for (j in 1:length(row_dims)) {
            data_set[[j]] <- abs(data$data_views[[j]]) + abs(mvrnorm(
                n = row_dims[j],
                mu = rep(0, col_dims[j]), Sigma = diag(noise_vec[i], col_dims[j])
            ))
        }
        openxlsx::write.xlsx(data_set, file = paste0(file_path2, "/data.xlsx")) # nolint
        openxlsx::write.xlsx(data$truth_rows, file = paste0(file_path2, "/true_rows.xlsx")) # nolint: line_length_linter.
        openxlsx::write.xlsx(data$truth_cols, file = paste0(file_path2, "/true_cols.xlsx")) # nolint: line_length_linter.
    }
}
