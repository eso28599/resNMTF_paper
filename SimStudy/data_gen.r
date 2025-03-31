args <- commandArgs(trailingOnly = TRUE)
path_to_sim_folder <- as.character(args[1])
batch_folder <- as.character(args[2])
source(paste0(
  path_to_sim_folder,
  "/sim_parameters.r"
)) # parameters for simulation
source("Functions/data_generation.r") # data gen files

path_to_data_folder <- paste0(
  paste0(path_to_sim_folder, "/data/"),
  batch_folder
)

if (factor == "noise") {
  save_data_noise(row_cl_dims, col_cl_dims, k_vals[1], path_to_data_folder,
    method_vec, noise_vec,
    row_same_shuffle = row_same_shuffle, col_same_shuffle = col_same_shuffle
  )
} else if (factor == "signal") {
  for (i in seq_len(length(method_vec))) {
    # generate data and save in given path
    path <- paste0(path_to_data_folder, paste0("/", method_vec[i]))
    save_data(row_cl_dims[i][[1]], col_cl_dims[i][[1]],
      k_vals[i],
      path, noise_level,
      row_same_shuffle = row_same_shuffle,
      col_same_shuffle = col_same_shuffle,
      signal = signal_vec[i]
    )
  }
} else if (factor == "overlap") {
  for (i in seq_len(length(method_vec))) {
    # generate data and save in given path
    path <- paste0(path_to_data_folder, paste0("/", method_vec[i]))
    save_data(row_cl_dims[i][[1]], col_cl_dims[i][[1]],
      k_vals[i],
      path, noise_level, row_e[i], col_e[i], row_o[i], col_o[i],
      row_same_shuffle = row_same_shuffle, col_same_shuffle = col_same_shuffle
    )
  }
} else {
  for (i in seq_len(length(method_vec))) {
    # generate data and save in given path
    path <- paste0(path_to_data_folder, paste0("/", method_vec[i]))
    save_data(row_cl_dims[i][[1]], col_cl_dims[i][[1]],
      k_vals[i],
      path, noise_level,
      row_same_shuffle = row_same_shuffle, col_same_shuffle = col_same_shuffle
    )
  }
}
