source("Functions/extra_funcs.r")
#file is specific to simulation
#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
n_views <- 3
row_cl_dims <- rep(200, n_views)
#100, 50,250 features respectively
col_cl_dims <- c(100, 50, 250)
# col_cl_dims <- c(200, 150, 250)
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

niter <- 200
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3)] <- 1
phi[2, c(3)] <- 1
val <- 200
phi_mat <- val * phi
overlap <- seq(0,40, 5)
row_o <- overlap/100
col_o <- overlap/100
row_e <- rep(1,length(overlap))
col_e <- rep(1,length(overlap))
k_vals <- rep(5, length(overlap))
method_vec <- paste0("/res_nmtf_", overlap)
row_cl_dims <- rep(list(row_cl_dims), length(method_vec))
col_cl_dims <- rep(list(col_cl_dims), length(method_vec))
method_vec_sgl <- paste0("/nmtf_", overlap)
method_vec_gfa <- paste0("/gfa_", overlap)
method_vec_issvd <- paste0("/issvd_", overlap)
factor <- "overlap"
kept_factor <- paste0(n_views, " views, ") #changed
factor_name <- "overlap"
x_title <- "Overlap rate"
plot_title <- "The effect of increasing the number of biclusters on performance"
file_names <- c(paste0("res_nmtf_", overlap),
    paste0("gfa_", overlap), 
    paste0("issvd_", overlap), paste0("nmtf_", overlap))
n_vars <- length(overlap)
method_vec_res <- c(rep("ResNMTF", n_vars), rep("GFA", n_vars),
 rep("iSSVD", n_vars), rep("NMTF", n_vars))
# method_vec_res <- c(rep("ResNMTF", 4),rep("NMTF", 4))
factor_vec <- rep(overlap, 4)
k_vec <- rep(overlap, 4)
phi_constant <- TRUE
n_col_plots <- 2
n_row_plots <- 2
col_names_tables <- overlap
order_fac <- FALSE 
