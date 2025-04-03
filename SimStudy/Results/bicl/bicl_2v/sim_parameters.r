source("Functions/extra_funcs.r")
# file is specific to simulation
# simulation parameters
noise_level <- 5
# each dataset with 200 samples, same clusters across views
n_views <- 2
row_cl_dims <- rep(200, n_views)
# 100, 50,250 features respectively
col_cl_dims <- c(100, 250)
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

k_vals <- c(3, 4, 5, 6)

niter <- 200
# parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, 2] <- 1
val <- 200
phi_mat <- val * phi
bicl_numbers <- 3:6
method_vec <- c(
  "/res_nmtf_3",
  "/res_nmtf_4", "/res_nmtf_5", "/res_nmtf_6"
)
row_cl_dims <- rep(list(row_cl_dims), length(method_vec))
col_cl_dims <- rep(list(col_cl_dims), length(method_vec))
method_vec_sgl <- paste0("/nmtf_", bicl_numbers)
method_vec_gfa <- paste0("/gfa_", bicl_numbers)
method_vec_issvd <- paste0("/issvd_", bicl_numbers)
factor <- "bicl"
kept_factor <- paste0(n_views, " views, ") # changed
factor_name <- "Biclusters"
x_title <- "Number of biclusters"
plot_title <- "The effect of increasing the number of biclusters on performance"
file_names <- c(
  paste0("res_nmtf_", bicl_numbers),
  paste0("gfa_", bicl_numbers),
  paste0("issvd_", bicl_numbers), paste0("nmtf_", bicl_numbers)
)
method_vec_res <- c(
  rep("ResNMTF", 4), rep("GFA", 4),
  rep("iSSVD", 4), rep("NMTF", 4)
)
factor_vec <- rep(3:6, 4)
k_vec <- rep(3:6, 4)
phi_constant <- TRUE
n_col_plots <- 2
n_row_plots <- 2
col_names_tables <- 3:6
order_fac <- FALSE
