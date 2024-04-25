args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1]) #simulation folder
batch_folder = as.character(args[2]) #repeat
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
source("resNMTF_funcs.r")

data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
data_views <- import_matrix(paste0(data_name, "/data.xlsx"), col_n = TRUE)
for (i in 1:length(method_vec)){
    #read in data
    file_path <- paste0(data_name, method_vec[i])
    #data_views <- import_matrix(paste0(file_path, "/data.xlsx"), col_n = TRUE)
    #nmtf
    nmtf_results <- restMultiNMTF_run(Xinput = data_views, phi = phi_mat, stab_thres = stab_vec[i])
    save_results(nmtf_results, file_path)
}
