source("Functions/extra_funcs.r")
#file is specific to simulation
#simulation parameters
noise_level <- 10
#each dataset with 200 samples, same clusters across views
n_views <- 5
row_cl_dims <- list(rep(200, 2), rep(200, 3),
                rep(200, 4), rep(200, 5))
#100, 50,250 features respectively
col_cl_dims <- list(rep(150, 2), rep(150, 3),
                rep(150, 4), rep(150, 5))
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

k_vals <- rep(5, 5)

niter <- 200
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3, 4, 5)] <- 1
phi[2, c(3, 4, 5)] <- 1
phi[3, c(4, 5)] <- 1
phi[4, 5] <- 1
val <- 200
phi_mat <- val * phi
phi_mat <- list(phi_mat[1:2,1:2], phi_mat[1:3,1:3], phi_mat[1:4,1:4], phi_mat)
bicl_numbers <- 2:5
method_vec <- c("/res_nmtf_2", "/res_nmtf_3",
 "/res_nmtf_4", "/res_nmtf_5")
# row_cl_dims <- rep(list(row_cl_dims), length(method_vec))
# col_cl_dims <- rep(list(col_cl_dims), length(method_vec))
method_vec_sgl <- paste0("/nmtf_", bicl_numbers)
method_vec_gfa <- paste0("/gfa_", bicl_numbers)
method_vec_issvd <- paste0("/issvd_", bicl_numbers)
factor <- "views"
kept_factor <- paste0(n_views, " views, ") #changed
factor_name <- "Views"
x_title <- "Number of views"
plot_title <- "The effect of increasing the number of views on performance"
file_names <- c(paste0("res_nmtf_", bicl_numbers),
    paste0("gfa_", bicl_numbers), 
    paste0("issvd_", bicl_numbers), paste0("nmtf_", bicl_numbers))
method_vec_res <- c(rep("ResNMTF", 4), rep("GFA", 4),
 rep("iSSVD", 4), rep("NMTF", 4))
factor_vec <- rep(2:5, 4)
k_vec <- rep(2:5, 4)
phi_constant <- TRUE
n_col_plots <- 2
n_row_plots <- 2
col_names_tables <- 2:5
order_fac <- FALSE 
