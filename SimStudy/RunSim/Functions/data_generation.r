library(openxlsx)
library(MASS)
library(Matrix)
##### synthetic data generation

#preparation of linked dataviews
# function to generate a dataview for specified bicluster dimensions
# same number of row and column clusters 
# row clusters remain common across all views, columns change 
get_dims <- function(n_r, n_c, k){
    part_r <- floor(n_r/7)
    part_c <- floor(n_c/7)
    if(k == 2){
        rows <- c(4 * part_r, n_r - 4 * part_r)
        cols <- c(4 * part_c, n_c - 4 * part_c)
    }else if (k == 3) {
        rows <- c(3 * part_r, 2 * part_r, n_r - 5 * part_r)
        cols <- c(3 * part_c, 2 * part_c, n_c - 5 * part_c)
    }else if (k == 4) {
        rows <- c(3 * part_r, 2 * part_r, part_r,  n_r - 6 * part_r)
        cols <- c(3 * part_c, 2 * part_c, part_c,  n_c - 6 * part_c)
    }else if (k == 5) {
        rows <- c(2 * part_r, 2 * part_r, part_r, part_r,  n_r - 6 * part_r)
        cols <- c(2 * part_c, 2 * part_c, part_c, part_c,  n_c - 6 * part_c)
    }else if (k == 6) {
        rows <- c(2 * part_r, part_r, part_r, part_r, part_r,  n_r - 6 * part_r)
        cols <- c(2 * part_c, part_c, part_c, part_c, part_c,  n_c - 6 * part_c)
    }
    rows <- sort(rows, decreasing = TRUE)
    cols <- sort(cols, decreasing = TRUE)
    row_start <- c(1, cumsum(rows)[-k] + 1)
    row_end <- cumsum(rows)
    col_start <- c(1, cumsum(cols)[-k] + 1)
    col_end <- cumsum(cols)
    
    return(list("Row_s" = row_start, "Row_e" = row_end,
             "Col_s" = col_start, "Col_e" = col_end))
}

one_view_adv <- function(row_dims, col_dims, k, noise,signal, overlap = FALSE){
    #'row_dims: vector of sizes of each row cluster in this view
    #'col_dims: vector of sizes of each row cluster in this view
    #'noise: variance of the noise added to views
    #' 
    #'view: matrix of this dataview
    #'truth_row: vector indicating membership of row clusters
    #'truth_col: vector indicating membership of col clusters

    # generate noise
    X_noise <- mvrnorm(n = row_dims, mu = rep(0, col_dims), Sigma = diag(noise, col_dims))

    # Create a list as input for block-diagonal data generation
    # list length of no of clusters in this first view

    # generate a mvn with n=number of individuals of view and no of features
    # equal to number of features of the view
    # mean of 5 for each column
    # covariance matrix is identity - each feature is independent
    #define vectors indicating true row/column cluster membership for each view
    true_row <- matrix(0, nrow = row_dims, ncol = k)
    true_col <- matrix(0, nrow = col_dims, ncol = k)
    X_view <- matrix(0, nrow = row_dims, ncol = col_dims)
    dims <- get_dims(row_dims, col_dims, k)
    row_start <- dims$Row_s
    row_end <- dims$Row_e
    col_start <- dims$Col_s
    col_end <- dims$Col_e
    #row_start
    for (i in 1:k){
        n_r <- (row_end[i] - row_start[i] + 1)
        n_c <- (col_end[i] - col_start[i] + 1)
        X_view[(row_start[i]):(row_end[i]), (col_start[i]):(col_end[i])] <-  mvrnorm(n = n_r, mu = rep(signal, n_c), Sigma = diag(n_c))
        # X_view[(row_start[i]):(row_end[i]), (col_start[i]):(col_end[i])] <-  mvrnorm(n = n_r, mu = rep(100, n_c), Sigma = diag(n_c))

        true_row[(row_start[i]):(row_end[i]), i] <- 1
        true_col[(col_start[i]):(col_end[i]), i] <- 1
    }
    # add noise to first view
    X_view <- abs(X_view) + abs(X_noise)
    return(list(view = X_view, truth_row = true_row, truth_col = true_col))
}

multi_view <- function(row_dims, col_dims, k, noise, signal, row_same_shuffle=TRUE, col_same_shuffle=TRUE, seed=FALSE){
    #'rowClusters: n length list of vectors of row cluster sizes in each view
    #'rowClusters: n length list of vectors of column cluster sizes in each view
    #'seed: logical indicator, default is FALSE, if true sets seed so same data is generated each time
    #' 
    #'data_views: n length list of data views
    #'truth_rows: n length list of vectors of true row cluster membership for each view
    #'col_rows: n length list of vectors of true column cluster membership for each view
    
    if(is.numeric(seed)){
        set.seed(seed)
    }
    n_views <- length(row_dims)
    # Introduce simulated data-views- where we store the views
    X_trial <- vector("list", length = n_views)
    #list to store row clusterings for each dataview
    true_row_clusterings <- vector("list", length = n_views)
    true_col_clusterings <- vector("list", length = n_views) 
    #initialise index to shuffle rows and columns 
    #new row/col index
    #if shuffle is true - dims of each view must be the same and represent same objects
    if(row_same_shuffle){
        new_row_ind <- sample(row_dims[1])
    }
    if(col_same_shuffle){
        new_col_ind <- sample(col_dims[1])
    }
    #generate data and store
    for(i in 1:n_views){
        if(!row_same_shuffle){
            new_row_ind <- sample(row_dims[i])
        }
        if(!col_same_shuffle){
            new_col_ind <- sample(col_dims[i])
        }
        data_i <- one_view_adv(row_dims[i], col_dims[i],
                k, noise, signal)
        X_trial[[i]] <- (data_i$view)[new_row_ind, new_col_ind]
        true_row_clusterings[[i]] <- (data_i$truth_row)[new_row_ind, ]
        true_col_clusterings[[i]] <- (data_i$truth_col)[new_col_ind, ]
    }

    return(list(data_views= X_trial, truth_rows = true_row_clusterings, truth_cols = true_col_clusterings))
}

save_data <- function(row_dims, col_dims, k, file_path, noise, row_same_shuffle=TRUE,col_same_shuffle=TRUE, signal=5){
        #can change the noise parameter here for level of noise in views
        data <- multi_view(row_dims, col_dims,k, noise, signal, row_same_shuffle , col_same_shuffle)
        #save data as a file in given directory
        #export each data frame to separate sheets in same Excel file
        openxlsx::write.xlsx(data$data_views, file = paste0(file_path, "/data.xlsx")) # nolint
        openxlsx::write.xlsx(data$truth_rows, file = paste0(file_path, "/true_rows.xlsx")) # nolint: line_length_linter.
        openxlsx::write.xlsx(data$truth_cols, file = paste0(file_path, "/true_cols.xlsx")) # nolint: line_length_linter.
}


save_data_noise <- function(row_dims, col_dims, k, file_path, method_vec, noise_vec, row_same_shuffle=TRUE,col_same_shuffle=TRUE){
        #can change the noise parameter here for level of noise in views
        data <- multi_view(row_dims, col_dims, k, 0, row_same_shuffle , col_same_shuffle )
        #save data as a file in given directory
        #export each data frame to separate sheets in same Excel file
        for(i in 1:length(noise_vec)){
            file_path2 <- paste0(file_path, paste0("/", method_vec[i]))
            data_set <- vector("list", length=length(row_dims))
            for(j in 1:length(row_dims)){
                data_set[[j]] <- abs(data$data_views[[j]]) + abs(mvrnorm(n = row_dims[j],
                             mu = rep(0, col_dims[j]), Sigma = diag(noise_vec[i], col_dims[j])))
            }   
            openxlsx::write.xlsx(data_set, file = paste0(file_path2, "/data.xlsx")) # nolint
            openxlsx::write.xlsx(data$truth_rows, file = paste0(file_path2, "/true_rows.xlsx")) # nolint: line_length_linter.
            openxlsx::write.xlsx(data$truth_cols, file = paste0(file_path2, "/true_cols.xlsx")) # nolint: line_length_linter.
        }     
}
