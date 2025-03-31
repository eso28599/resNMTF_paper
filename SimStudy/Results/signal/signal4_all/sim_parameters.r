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

k_vals <- rep(4, 20)
signal_vec <- (1:20)*5
niter <- 200
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3)] <- 1
phi[2, c(3)] <- 1
val <- 200
phi_mat <- val * phi
bicl_numbers <- 5
param_vals <- 1:20
method_vec <- paste0("/res_nmtf_", param_vals)
row_cl_dims <- rep(list(row_cl_dims), length(method_vec))
col_cl_dims <- rep(list(col_cl_dims), length(method_vec))
method_vec_issvd <- paste0("/issvd_", param_vals)
method_vec_sgl <- paste0("/nmtf_", param_vals)
method_vec_gfa <- paste0("/gfa_", param_vals)
method_vec_issvd <- paste0("/issvd_", param_vals)
factor <- "signal"
kept_factor <- paste0(n_views, " views, ") #changed
factor_name <- "Parameter"
x_title <- "Signal"
plot_title <- "The effect of increasing the signal of biclusters on performance"
file_names <- c(paste0("res_nmtf_", param_vals),
    paste0("gfa_", param_vals), 
    paste0("issvd_", param_vals), paste0("nmtf_", param_vals))

method_vec_res <- c(rep("ResNMTF", 20), rep("GFA", 20),
 rep("iSSVD", 20), rep("NMTF", 20))
factor_vec <- rep((1:20)*5, 4)
k_vec <-rep(5,80)
phi_constant <- TRUE
n_col_plots <- 2
n_row_plots <- 2
col_names_tables <- (1:20)*5
order_fac <- FALSE 
