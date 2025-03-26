args <- commandArgs(trailingOnly = TRUE)
path_to_sim_folder <- as.character(args[1])
batch_folder <- as.character(args[2])

# simulation evaluation file
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
source("Functions/evaluation_funcs.r")

if (factor == "signal") {
  # analyse each issvd parameter result
  for (i in seq_along(method_vec_issvd)) {
    # read in results
    data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
    true_rows <- import_matrix(
      paste0(data_name, "/issvd_", i, "/true_rows.xlsx")
    )
    true_cols <- import_matrix(
      paste0(data_name, "/issvd_", i, "/true_cols.xlsx")
    )
    data_views <- import_matrix(
      paste0(data_name, "/issvd_", i, "/data.xlsx")
    )
    eval_method(
      data_name, method_vec_issvd[i],
      true_rows, true_cols, data_views
    )
  }
} else {
  # read in results
  data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
  true_rows <- import_matrix(paste0(data_name, "/true_rows.xlsx"))
  true_cols <- import_matrix(paste0(data_name, "/true_cols.xlsx"))
  data_views <- import_matrix(paste0(data_name, "/data.xlsx"))
  # analyse each result
  for (i in seq_along(method_vec_issvd)) {
    eval_method(
      data_name, method_vec_issvd[i],
      true_rows, true_cols, data_views
    )
  }
}
