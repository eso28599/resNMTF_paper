source("Functions/extra_funcs.r")
#file is specific to simulation
#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
n_views <- 2
row_cl_dims <- rep(100, n_views)
#100, 50,250 features respectively
col_cl_dims <- c(1000, 1000)
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

k_vals <- 4

niter <- 200
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, 2] <- 1
val <- 200
phi_mat <- val * phi
bicl_numbers <- 5
param_vals <- 1:10
sigma_vec <- (1:10)/10
# method_vec <- paste0("/issvd_", param_vals)
# method_vec <- paste0("")
method_vec <- paste0("/res_nmtf_", param_vals)
row_cl_dims <- rep(list(row_cl_dims), length(method_vec))
col_cl_dims <- rep(list(col_cl_dims), length(method_vec))
# row_cl_dims <- rep(list(row_cl_dims), length(method_vec))
# col_cl_dims <- rep(list(col_cl_dims), length(method_vec))
method_vec_sgl <- paste0("/nmtf_", param_vals)
method_vec_gfa <- paste0("/gfa_", param_vals)
method_vec_issvd <- paste0("/issvd_", param_vals)
factor <- "issvd_param"
kept_factor <- paste0(n_views, " views, ") #changed
factor_name <- "Parameter"
x_title <- "sigma"
plot_title <- "The effect of increasing sigma on performance"
file_names <- paste0("issvd_", param_vals)
method_vec_res <- rep("iSSVD", 10)
sigma_vec <- (1:10)/10
k_vec <- rep(4,10)
phi_constant <- FALSE
n_col_plots <- 2
n_row_plots <- 2
col_names_tables <- 1:10
order_fac <- FALSE 
