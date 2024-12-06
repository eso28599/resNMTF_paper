source("Functions/extra_funcs.r")
#file is specific to simulation
#simulation parameters
noise_vec <- c(1:9,10*(1:10))
#each dataset with 200 samples, same clusters across views
n_views <- 3
row_cl_dims <-  rep(200, n_views)
#100, 50,250 features respectively
col_cl_dims <- c(100, 50, 250)
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

k_vals <- rep(5, 19)

niter <- 200
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3)] <- 1
phi[2, c(3)] <- 1
val <- 200
phi_mat <- val * phi
bicl_numbers <- c(1:9,10*(1:10))
method_vec <- paste0("/res_nmtf_", c(1:9,10*(1:10)))
row_cl_dims <- rep(list(row_cl_dims), length(method_vec))
col_cl_dims <- rep(list(col_cl_dims), length(method_vec))
method_vec_sgl <- paste0("/nmtf_", bicl_numbers)
method_vec_gfa <- paste0("/gfa_", bicl_numbers)
method_vec_issvd <- paste0("/issvd_", bicl_numbers)
factor <- "noise"
kept_factor <- paste0(n_views, " views, ") #changed
factor_name <- "Noise"
x_title <- "Level of noise"
plot_title <- "The effect of increasing the level of noise on performance"
file_names <- c(paste0("res_nmtf_", bicl_numbers),
    paste0("gfa_", bicl_numbers), 
    paste0("issvd_", bicl_numbers), paste0("nmtf_", bicl_numbers))
method_vec_res <- c(rep("ResNMTF", 19), rep("GFA", 19),
 rep("iSSVD", 19), rep("NMTF", 19))
factor_vec <- rep(c(1:9,10*(1:10)), 4)
k_vec <- rep(c(1:9,10*(1:10)), 4)
phi_constant <- TRUE
n_col_plots <- 2
n_row_plots <- 2
col_names_tables <- c(1:9,10*(1:10))
order_fac <- FALSE 
