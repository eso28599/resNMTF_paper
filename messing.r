
# checking the results of the issvd on data generated from it's approach
ssvd_data2 <- list(X1 + mvrnorm(100, rep(0, 1000), 0.01*diag(1, 1000)),
        X2 + mvrnorm(100, rep(0, 1000), 0.01*diag(1, 1000)))
file_path <- "test_data"
openxlsx::write.xlsx(ssvd_data2, file = paste0(file_path, "/data.xlsx")) # nolint
openxlsx::write.xlsx(true_rows, file = paste0(file_path, "/true_rows.xlsx")) # nolint: line_length_linter.
openxlsx::write.xlsx(true_cols, file = paste0(file_path, "/true_cols.xlsx")) #

#check issvd results 
source("SimStudy/RunSim/Functions/extra_funcs.r")
source("SimStudy/RunSim/Functions/evaluation_funcs.r")
source("SimStudy/RunSim/main.r")
true_rows <- import_matrix("test_data/true_rows.xlsx")
true_cols <- import_matrix("test_data/true_cols.xlsx")
data <- import_matrix("test_data/data.xlsx")
rows_issvd <- import_matrix("test_data/issvd_rows.xlsx")
cols_issvd <- import_matrix("test_data/issvd_cols.xlsx")
jaccard_res(rows_issvd[[2]], cols_issvd[[2]], true_rows[[2]], true_cols[[2]])
jaccard_res(rows_issvd[[1]], cols_issvd[[1]], true_rows[[1]], true_cols[[1]])

# applying the restMultiNMTF on the data generated from issvd
phi_mat <- matrix(0,2,2)
phi_mat[1,2] <- 200
res_og2 <- restMultiNMTF_run(data, phi=2*phi_mat)
jaccard_res(res_og2$row_clusters[[1]], res_og2$col_clusters[[1]], true_rows[[1]], true_cols[[1]])
jaccard_res(res_og2$row_clusters[[2]], res_og2$col_clusters[[2]], true_rows[[2]], true_cols[[2]])


# generate synthetic data via our approach. 
source("main.r")
source('SimStudy/RunSim/Functions/data_generation.r')
# generate a test dataset
n_views <- 5
row_cl_dims <- rep(200, n_views)
#100, 50,250 features respectively
col_cl_dims <- c(100, 50, 250)
col_cl_dims <- rep(250,5)
save_data(row_cl_dims, col_cl_dims, 5, 'test_data2', 5 ,col_same_shuffle=FALSE, signal=5)
data <- import_matrix("test_data2/data.xlsx")
rows <- import_matrix("test_data2/true_rows.xlsx")
cols <- import_matrix("test_data2/true_cols.xlsx")


niter <- 200
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3, 4, 5)] <- 1
phi[2, c(3, 4, 5)] <- 1
phi[3, c(4, 5)] <- 1
phi[4, 5] <- 1
val <- 200
res_og <- restMultiNMTF_run(data, phi=val*phi)
jaccard_res(res_og$row_clusters[[1]], res_og$col_clusters[[1]], rows[[1]], cols[[1]])
jaccard_res(res_og$row_clusters[[2]], res_og$col_clusters[[2]], rows[[2]], cols[[2]])

max(data[[1]])
#perform gfa
source('SimStudy/RunSim/OtherMethods/gfa_funcs.r')
source("SimStudy/RunSim/Functions/extra_funcs.r")
test <- gfa_apply(data, 10)
test_high <- gfa_apply(data, 10, 0.9)
# check results
true_rows <- import_matrix("test_data2/true_rows.xlsx")
true_cols <- import_matrix("test_data2/true_cols.xlsx")
jaccard_res(test_high$row_clusters[[1]], test_high$col_clusters[[1]], true_rows[[1]], true_cols[[1]])
jaccard_res(test_high$row_clusters[[2]], test_high$col_clusters[[2]], true_rows[[2]], true_cols[[2]])
jaccard_res(test$row_clusters[[1]], test$col_clusters[[1]], true_rows[[1]], true_cols[[1]])
jaccard_res(test$row_clusters[[2]], test$col_clusters[[2]], true_rows[[2]], true_cols[[2]])


rows_issvd <- import_matrix("test_data2/issvd_rows.xlsx")
cols_issvd <- import_matrix("test_data2/issvd_cols.xlsx")
jaccard_res(rows_issvd[[2]], cols_issvd[[2]], true_rows[[2]], true_cols[[2]])
jaccard_res(rows_issvd[[1]], cols_issvd[[1]], true_rows[[1]], true_cols[[1]])

# X2 <- (U %*% S %*% t(V2))
#reorder rows
ordered_rows <- c(which(true_rows[[1]][,1] == 1),
    which(true_rows[[1]][,2] == 1),
    which(true_rows[[1]][,3] == 1),
    which(true_rows[[1]][,4] == 1),
    which(rowSums(true_rows[[1]])==0))

ordered_cols <- c(which(true_cols[[1]][,1] == 1),
    which(true_cols[[1]][,2] == 1),
    which(true_cols[[1]][,3] == 1),
    which(true_cols[[1]][,4] == 1),
    which(rowSums(true_cols[[1]])==0))

ordered_cols2 <- c(which(true_cols[[2]][,1] == 1),
    which(true_cols[[2]][,2] == 1),
    which(true_cols[[2]][,3] == 1),
    which(true_cols[[2]][,4] == 1),
    which(rowSums(true_cols[[2]])==0))

image(t(X1[ordered_rows, ordered_cols]+ mvrnorm(100, rep(0, 1000), 1*diag(1, 1000))))
image(t(X2[ordered_rows, ordered_cols2]))

heatmap(t(X1+ mvrnorm(100, rep(0, 1000), 1*diag(1, 1000))))

path_to_save <- paste0(paste0(path_to_sim_folder, "/data/"),
                         batch_folder)

ssvd_data <- list(X1 + mvrnorm(100, rep(0, 1000), 1*diag(1, 1000)),
        X2 + mvrnorm(100, rep(0, 1000), 1*diag(1, 1000)))




source("main.r")
res_og <- restMultiNMTF_run(ssvd_data, phi=phi_mat, stability = FALSE)
res_og2 <- restMultiNMTF_run(ssvd_data2, phi=phi_mat)


source("SimStudy/Functions/evaluation_funcs.r")
source("SimStudy/Functions/extra_funcs.r")

jaccard_res(res_og$row_clusters[[1]], res_og$col_clusters[[1]], true_rows[[1]], true_cols[[1]])
jaccard_res(res_og$row_clusters[[1]], res_og$col_clusters[[1]], true_rows[[1]], true_cols[[1]])
jaccard_res(res_og2$row_clusters[[1]], res_og2$col_clusters[[1]], true_rows[[1]], true_cols[[1]])




rows_issvd <- import_matrix("test_data/issvd_rows.xlsx")
cols_issvd <- import_matrix("test_data/issvd_cols.xlsx")
true_rows <- import_matrix("test_data/true_rows.xlsx")
true_cols <- import_matrix("test_data/true_cols.xlsx")



