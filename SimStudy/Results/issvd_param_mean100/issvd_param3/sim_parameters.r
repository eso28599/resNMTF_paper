source("Functions/extra_funcs.r")
#file is specific to simulation
#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
n_views <- 3
row_cl_dims <- rep(200, n_views)
#100, 50,250 features respectively
col_cl_dims <- c(100, 50, 250)
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

k_vals <- 3

niter <- 200
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3)] <- 1
phi[2, c(3)] <- 1
val <- 200
phi_mat <- val * phi
bicl_numbers <- 5
param_vals <- 0:20
# method_vec <- paste0("/issvd_", param_vals)
method_vec <- paste0("")
row_cl_dims <- rep(list(row_cl_dims), length(method_vec))
col_cl_dims <- rep(list(col_cl_dims), length(method_vec))
method_vec_issvd <- paste0("/issvd_", param_vals)
factor <- "issvd_param"
kept_factor <- paste0(n_views, " views, ") #changed
factor_name <- "Parameter"
x_title <- "Per-comparrison wise error rate"
plot_title <- "The effect of increasing the per-comparrison wise error rate on iSSVD performance"
file_names <- paste0("issvd_", param_vals)
method_vec_res <- rep("iSSVD", 21)
factor_vec <- 0:21/20
k_vec <- 0:21
phi_constant <- FALSE
n_col_plots <- 2
n_row_plots <- 2
col_names_tables <- 3:6
order_fac <- FALSE 
