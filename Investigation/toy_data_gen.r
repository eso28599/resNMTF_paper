library(openxlsx)
library(MASS)
library(Matrix)

#toy examples
get_dims_toy <- function(n_r, n_c, k, row_e = 1, col_e = 1, row_o = 0, col_o = 0){
    #' row_e portion of rows in biclusters
    part_r <- floor((row_e*n_r)/7)
    part_c <- floor((col_e*n_c)/7)
    if(row_e == 1){
        if(k == 2){
            rows <- c(4 * part_r, n_r - 4 * part_r)
        }else if (k == 3) {
            rows <- c(3 * part_r, 2 * part_r, n_r - 5 * part_r)
        }else if (k == 4) {
            rows <- c(3 * part_r, 2 * part_r, part_r,  n_r - 6 * part_r)
        }else if (k == 5) {
            rows <- c(2 * part_r, 2 * part_r, part_r, part_r,  n_r - 6 * part_r)
        }else if (k == 6) {
            rows <- c(2 * part_r, part_r, part_r, part_r, part_r,  n_r - 6 * part_r)
        }
    }else{
       if(k == 2){
            rows <- c(4 * part_r, 3 * part_r)
        }else if (k == 3) {
            rows <- c(3 * part_r, 2 * part_r, 2 * part_r)
        }else if (k == 4) {
            rows <- c(3 * part_r, 2 * part_r, part_r,  part_r)
        }else if (k == 5) {
            rows <- c(2 * part_r, 2 * part_r, part_r, part_r,  part_r)
        }else if (k == 6) {
            rows <- c(2 * part_r, part_r, part_r, part_r, part_r,  part_r)
        } 
    }

    if(col_e == 1){
        if(k == 2){
            cols <- c(4 * part_c, n_c - 4 * part_c)
        }else if (k == 3) {
            cols <- c(3 * part_c, 2 * part_c, n_c - 5 * part_c)
        }else if (k == 4) {
            cols <- c(3 * part_c, 2 * part_c, part_c,  n_c - 6 * part_c)
        }else if (k == 5) {
            cols <- c(2 * part_c, 2 * part_c, part_c, part_c,  n_c - 6 * part_c)
        }else if (k == 6) {
            cols <- c(2 * part_c, part_c, part_c, part_c, part_c,  n_c - 6 * part_c)
        }
    }else{
       if(k == 2){
            cols <- c(4 * part_c, 3 * part_c)
        }else if (k == 3) {
            cols <- c(3 * part_c, 2 * part_c, 2 * part_c)
        }else if (k == 4) {
            cols <- c(3 * part_c, 2 * part_c, part_c,  part_c)
        }else if (k == 5) {
            cols <- c(2 * part_c, 2 * part_c, part_c, part_c,  part_c)
        }else if (k == 6) {
            cols <- c(2 * part_c, part_c, part_c, part_c, part_c,  part_c)
        } 
    }

    
    rows <- sort(rows, decreasing = TRUE)
    cols <- sort(cols, decreasing = TRUE)

    row_overlap <- floor(rows*row_o)
    row_overlap[k] <- 0
    col_overlap <- floor(cols*col_o)
    col_overlap[k] <- 0

    row_start <- c(1, cumsum(rows)[-k] + 1)
    if(part_r==0){
        row_start <- rep(0, k)
    }
    row_end <- cumsum(rows)
    row_end <- row_end + row_overlap
    col_start <- c(1, cumsum(cols)[-k] + 1)
    if(part_c==0){
        col_start <- rep(0, k)
    }
    col_end <- cumsum(cols)
    col_end <- col_end + col_overlap
    
    return(list("Row_s" = row_start, "Row_e" = row_end,
             "Col_s" = col_start, "Col_e" = col_end))
}



one_view_toy <- function(row_dims, col_dims, k, noise, row_e = 1, col_e = 1, row_o = 0, col_o = 0){
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
    dims <- get_dims_toy(row_dims, col_dims, k, row_e, col_e, row_o, col_o)
    row_start <- dims$Row_s
    row_end <- dims$Row_e
    col_start <- dims$Col_s
    col_end <- dims$Col_e
    #row_start
    for (i in 1:k){
        n_r <- (row_end[i] - row_start[i] + 1)
        n_c <- (col_end[i] - col_start[i] + 1)
        if(n_r !=0){
            X_view[(row_start[i]):(row_end[i]), (col_start[i]):(col_end[i])] <-  X_view[(row_start[i]):(row_end[i]), (col_start[i]):(col_end[i])] +
                    mvrnorm(n = n_r, mu = rep(5, n_c), Sigma = diag(n_c))
            true_row[(row_start[i]):(row_end[i]), i] <- 1
            true_col[(col_start[i]):(col_end[i]), i] <- 1
        }  
    }
    # add noise to first view
    X_view <- abs(X_view) + abs(X_noise)
    return(list(view = X_view, truth_row = true_row, truth_col = true_col))
}

make_longer <- function(vec, n){
    if(length(vec)==1){
        vec <- rep(vec, n)
    }
    return(vec)
}

multi_view_toy <- function(file_path, row_dims, col_dims, k_vec, noise_vec,row_e = 1, col_e = 1, row_o = 0 , col_o = 0, row_same_shuffle=TRUE, col_same_shuffle=FALSE, seed=FALSE){
    #'rowClusters: n length list of vectors of row cluster sizes in each view
    #'rowClusters: n length list of vectors of column cluster sizes in each view
    #'seed: logical indicator, default is FALSE, if true sets seed so same data is generated each time
    #' 
    #'data_views: n length list of data views
    #'truth_rows: n length list of vectors of true row cluster membership for each view
    #'col_rows: n length list of vectors of true column cluster membership for each view
    
    if(seed){
        set.seed(20)
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
    k_vec <- make_longer(k_vec, n_views)
    noise_vec <- make_longer(noise_vec, n_views)
    row_e <- make_longer(row_e, n_views)
    row_o <- make_longer(row_o, n_views)
    col_e <- make_longer(col_e, n_views)
    col_o <- make_longer(col_o, n_views)
    #generate data and store
    for(i in 1:n_views){
        if(!row_same_shuffle){
            new_row_ind <- sample(row_dims[i])
        }
        if(!col_same_shuffle){
            new_col_ind <- sample(col_dims[i])
        }
        data_i <- one_view_toy(row_dims[i], col_dims[i],
                k_vec[i], noise_vec[i], row_e[i], col_e[i], row_o[i], col_o[i])
        # image(t(data_i$view))
        X_trial[[i]] <- (data_i$view)[new_row_ind, new_col_ind]
        true_row_clusterings[[i]] <- (data_i$truth_row)[new_row_ind, ]
        true_col_clusterings[[i]] <- (data_i$truth_col)[new_col_ind, ]
        if(!is.na(file_path)){
            pdf(paste0(file_path, "/true_image_", i, ".pdf"))
            image(t(data_i$view))
            dev.off()
        }
    }
    #save data
    if(!is.na(file_path)){
        row_filename <- paste0(file_path, "/true_rows.xlsx")
        col_filename <- paste0(file_path, "/true_cols.xlsx")
        data_filename <- paste0(file_path, "/data.xlsx")
        openxlsx::write.xlsx(true_row_clusterings, file = row_filename)
        openxlsx::write.xlsx(true_col_clusterings, file = col_filename)
        openxlsx::write.xlsx(X_trial, file = data_filename)
    }
    
    return(list(data_views= X_trial, truth_rows = true_row_clusterings, truth_cols = true_col_clusterings))
}

choose_phi <- function(Xinput, phi_vals=NULL, Finput=NULL, Sinput=NULL, 
            Ginput=NULL, KK=NULL, phi=NULL, xi=NULL, psi=NULL, 
            nIter=NULL, k_min=3, k_max=6, repeats=5, no_clusts=FALSE, 
             sample_rate=0.8, n_stability=5, stability=TRUE, stab_thres=0.6){

                res_list <- vector("list", length=length(phi_vals))
                bisil_list <- rep(0, length=length(phi_vals))
                for(i in 1:length(phi_vals)){
                    res_list[[i]] <- restMultiNMTF_run(Xinput = Xinput,
                            Finput = Finput, 
                            Sinput = Sinput, 
                            Ginput = Ginput, 
                            KK = KK,
                            phi = phi_vals[i]*phi,
                            xi = xi, 
                            psi = psi, 
                            nIter = nIter, 
                            k_min = k_min, 
                            k_max = k_max,
                             repeats = repeats,
                            no_clusts = no_clusts, 
                            sample_rate = sample_rate, 
                            n_stability = n_stability, stability=stability,
                             stab_thres=stab_thres)
                    bisil_list[[i]] <- res_list[[i]]$Sil_score
                }
                return(list("scores" = bisil_list, "results" = res_list[[which.max(bisil_list)]], "all_res" = res_list))
}

toy_example <- function(data_obj, file_path, phi_vec=NULL, Finput=NULL, Sinput=NULL, 
            Ginput=NULL, KK=NULL, phi=NULL, xi=NULL, psi=NULL, 
            nIter=NULL, k_min=3, k_max=6, repeats=5, no_clusts=FALSE, 
             sample_rate=0.8, n_stability=5, stability=TRUE, stab_thres=0.6){
        nmtf <- choose_phi(Xinput = data_obj$data_views,
                            phi_vals=phi_vec,
                            phi = phi,
                            Finput = Finput, 
                            Sinput = Sinput, 
                            Ginput = Ginput, 
                            KK = KK,
                            xi = xi, 
                            psi = psi, 
                            nIter = nIter, 
                            k_min = k_min, 
                            k_max = k_max,
                             repeats = repeats,
                            no_clusts = no_clusts, 
                            sample_rate = sample_rate, 
                            n_stability = n_stability, stability=stability,
                             stab_thres=stab_thres)
    results <- vector("list", length=length(phi_vec))
    f_vec <- rep(0, length(phi_vec))
    rec_vec <- rep(0, length(phi_vec))
    rel_vec <- rep(0, length(phi_vec))
    for(i in 1:length(phi_vec)){
            results[[i]] <- evaluate_simulation_comp((nmtf$all_res[[i]])$row_clusters,
                 (nmtf$all_res[[i]])$col_clusters,
            data_obj$truth_rows, data_obj$truth_cols, data_obj$data_views)
            f_vec[i] <- (mean(results[[i]]$F_score))
            rec_vec[i] <- (mean(results[[i]]$Rec))
            rel_vec[i] <- (mean(results[[i]]$Rel))
            row_filename <- paste0(file_path, "/", phi_vec[i], "_row_clusts.xlsx")
            col_filename <- paste0(file_path, "/", phi_vec[i],"_col_clusts.xlsx")
            openxlsx::write.xlsx((nmtf$all_res[[i]])$row_clusters, file = row_filename)
            openxlsx::write.xlsx((nmtf$all_res[[i]])$col_clusters, file = col_filename)
    }
    max_b <- which.max(nmtf$scores)
    max_f <- which.max(f_vec)
    correct_choice <- max_f == max_b
    if(correct_choice){
        print(paste0("BiS correctly selects phi value corresponding to highest f score,", phi_vec[max_b],"."))
        print(paste0("This gives f-score ", max(f_vec),
         ", recovery score ", rec_vec[max_b],
          " and relevance score ", rel_vec[max_b],"."))
    }else{
        print(paste0("Selects ",
         phi_vec[max_b], " but ",phi_vec[max_f], " gives better performance." ))
        print(paste0("The selected phi value gives f-score ", f_vec[max_b],
         ", recovery score ", rec_vec[max_b],
          " and relevance score ", rel_vec[max_b],"."))
        print(paste0("The true optimal phi value gives f-score ", f_vec[max_f],
         ", recovery score ", rec_vec[max_f],
          " and relevance score ", rel_vec[max_f],"."))
    }
    summary <- data.frame("phi" = phi_vec, "f_score" = f_vec,
         "rec" =  rec_vec, "rel" =  rel_vec, "BiS" = nmtf$scores)
    write.csv(summary, file = paste0(file_path, "/summary.csv"))

    return(list("results" = results, "f_score" = f_vec, "rec" =  rec_vec, "rel" =  rel_vec, "BiS" = nmtf$scores))
}
