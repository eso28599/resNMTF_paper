args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
batch_folder = as.character(args[2])

# simulation evaluation file
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
source("Functions/evaluation_funcs.r")

#read in results
data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
#analyse each method
for (i in 1:length(method_vec)){
    file_path <- paste0(data_name, method_vec[i])
    #export each data frame to separate sheets in same Excel file
    true_rows <- import_matrix(paste0(file_path, "/true_rows.xlsx"))
    true_cols <- import_matrix(paste0(file_path, "/true_cols.xlsx"))
    data_views <- import_matrix(paste0(file_path, "/data.xlsx"))

    eval_method(data_name, method_vec[i])#resNMTF
    eval_method(data_name, method_vec_gfa[i])#GFA
    eval_method(data_name, method_vec_sgl[i])#NMTF
    eval_method(data_name, method_vec_issvd[i])#iSSVD
}
