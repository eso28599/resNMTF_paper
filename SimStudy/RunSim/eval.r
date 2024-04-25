args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
batch_folder = as.character(args[2])

# simulation evaluation file
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
source("evaluation_funcs.r")

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
    #now look at issvd
    results_rep <- vector(mode = "list", length = 5)
    sil_list <- c()
    file_path3 <- paste0(data_name, method_vec_issvd[i])
    for (j in 1:5){
        file <- paste0(file_path3, paste0("/", (j - 1)))
        row_filename <- paste0(file, "_row_clusts.xlsx")
        col_filename <- paste0(file, "_col_clusts.xlsx")
        #results for nmtf
        results_rep[[j]] <- evaluate_simulation_comp(import_matrix(row_filename),
                                                import_matrix(col_filename),
                                                true_rows, true_cols, data_views)
        sil_list <- c(sil_list, mean(results_rep[[j]]$BiS))
    }
    path_to_save <- paste0(file_path3, "_results.xlsx")
    openxlsx::write.xlsx(results_rep[[which.max(sil_list)]], file = path_to_save)
}
