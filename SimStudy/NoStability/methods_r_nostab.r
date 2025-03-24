args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1]) #simulation folder
batch_folder = as.character(args[2]) #repeat
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
source("OtherMethods/gfa_funcs.r")
source("main.r")

data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
for (i in 1:length(method_vec)){
    #read in data
    file_path <- paste0(data_name, method_vec[i])
    file_path2 <- paste0(data_name, method_vec_gfa[i])
    file_path3 <- paste0(data_name, method_vec_sgl[i])
    file_path4 <- paste0(data_name, method_vec_issvd[i])
    data_views <- import_matrix(paste0(file_path, "/data.xlsx"), col_n = TRUE)
    #nmtf
    if(factor == "views"){
        nmtf_results <- restMultiNMTF_run(Xinput = data_views, phi = phi_mat[i][[1]])
        save_results(nmtf_results, file_path)
        nmtf_sgl_res <- single_nmtf(data_views)
        save_results(nmtf_sgl_res, file_path3)
        # no stability
        nmtf_results0 <- restMultiNMTF_run(Xinput = data_views, phi = phi_mat[i][[1]], stability=FALSE)
        save_results(nmtf_results0, file_path2)
        nmtf_sgl_res0 <- single_nmtf(data_views, stab = FALSE)
        save_results(nmtf_sgl_res0, file_path4)
    }else{
        nmtf_results <- restMultiNMTF_run(Xinput = data_views, phi = phi_mat)
        nmtf_sgl_res <- single_nmtf(data_views)
    }   
}
