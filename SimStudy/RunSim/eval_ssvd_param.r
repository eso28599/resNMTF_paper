args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
batch_folder = as.character(args[2])

# simulation evaluation file
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
source("Functions/evaluation_funcs.r")

#read in results
data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
# for issvd param selection
true_rows <- import_matrix(paste0(data_name, "/true_rows.xlsx"))
true_cols <- import_matrix(paste0(data_name, "/true_cols.xlsx"))
data_views <- import_matrix(paste0(data_name, "/data.xlsx"))
#analyse each method
for (i in 1:length(method_vec_issvd)){
    eval_method(data_name, method_vec_issvd[i], true_rows, true_cols, data_views)#iSSVD

}
