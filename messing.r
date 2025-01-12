
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
source("main.r")

true_rows <- import_matrix("test_data/true_rows.xlsx")
true_cols <- import_matrix("test_data/true_cols.xlsx")
data <- import_matrix("test_data/data.xlsx")
rows_issvd <- import_matrix("test_data/issvd_rows.xlsx")
cols_issvd <- import_matrix("test_data/issvd_cols.xlsx")
jaccard_res(rows_issvd[[2]], cols_issvd[[2]], true_rows[[2]], true_cols[[2]])
jaccard_res(rows_issvd[[1]], cols_issvd[[1]], true_rows[[1]], true_cols[[1]])

# applying the restMultiNMTF on the data generated from issvd
phi_mat <- matrix(0,4,4)
phi_mat[1,c(2,3,4)] <- 200
phi_mat[2,c(3,4)] <- 200
phi_mat[3,c(4)] <- 200

phi_mat <- matrix(0,3,3)
phi_mat[1,c(2,3)] <- 200
phi_mat[2,3] <- 200


res_og_dg <- restMultiNMTF_run(data, phi=phi_mat, stability = FALSE, k_min=5, k_max=5)
res_og_stab <- restMultiNMTF_run(data, phi=phi_mat, stability = FALSE)
res_og2 <- restMultiNMTF_run(data, phi=phi_mat)
res_high_max <- restMultiNMTF_run(data, phi=phi_mat, k_max=8)

#old resNMTF function
res_old_main <- restMultiNMTF_run(data, phi=phi_mat, stab_thres = 0.4)
#old resNMTF function, with new clustering function
res_old_main_bis <- restMultiNMTF_run(data, phi=phi_mat, stab_thres = 0.4)
i<- 4
row_switch <- res_old_main_bis$Foutput[[i]]>1/200
col_switch <- res_old_main_bis$Goutput[[i]]>1/150
jaccard_res(row_switch, col_switch, true_rows[[i]], true_cols[[i]])


# problem in biclustering assignment/bisil functions!
results <- restMultiNMTF_run(data, phi=phi_mat,k_min=5, k_max=5, stability=FALSE)
bis <- bisil_score(data, results$row_clusters, results$col_clusters)
sil <- sil_score_inner(data[[1]], results$row_clusters[[1]], results$col_clusters[[1]])

plot(results$Foutput[[1]])

bis <- bisil_score(list(data[[1]]), list(row_switch), list(col_switch))

jaccard_res(row_switch, col_switch[,c(1,2,3,3)], true_rows[[i]], true_cols[[i]])


# 3v5b 200 features, unrestricted
res_0 <- restMultiNMTF_run(data, phi=0*phi_mat)
# 3v5b 200 features, restricted
res_200 <- restMultiNMTF_run(data, phi=phi_mat)


res_og_main <- restMultiNMTF_run(data, phi=phi_mat, stab_thres = 0.4)
res_og_main2 <- restMultiNMTF_run(data, phi=phi_mat, stab_thres = 0.4, sample_rate = 0.9)
res_og_main_oldbis <- restMultiNMTF_run(data, phi=phi_mat, stab_thres = 0.4)
res_og_main_oldbis <- restMultiNMTF_run(data, phi=phi_mat, stab_thres = 0.4)
res_og_oldmain <- restMultiNMTF_run(data, phi=phi_mat, stab_thres = 0.4)
res_og_just_oldmain <- restMultiNMTF_run(data, phi=phi_mat, stab_thres = 0.4)
res_og_oldmain_stats <- restMultiNMTF_run(data, phi=phi_mat, stab_thres = 0.4)
#new main, old get_thresholds, bisil_score, clustering_res_nmtf and main
res_og_old_main2 <- restMultiNMTF_run(data, phi=phi_mat, stab_thres = 0.4)

res_og40 <- restMultiNMTF_run(data, phi=phi_mat)
true_rows <- import_matrix("test_data2/true_rows.xlsx")
true_cols <- import_matrix("test_data2/true_cols.xlsx")
i<- 1
jaccard_res(res_og_stab$row_clusters[[i]], res_og_stab$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_og2$row_clusters[[i]], res_og2$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_200$row_clusters[[i]], res_200$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_og_main$row_clusters[[i]], res_og_main$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_og2b$row_clusters[[i]], res_og2b$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_og_main_oldbis$row_clusters[[i]], res_og_main_oldbis$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_og_oldmain$row_clusters[[i]], res_og_oldmain$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_og_just_oldmain$row_clusters[[i]], res_og_just_oldmain$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_og_oldmain_stats$row_clusters[[i]], res_og_oldmain_stats$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
accard_res(res_og2$Foutput[[2]]>1/200, res_og2$Goutput[[2]]>1/150, true_rows[[2]], true_cols[[2]])
jaccard_res(res_og_stab$Foutput[[i]]>1/200, res_og_stab$Goutput[[i]]>1/150, true_rows[[i]], true_cols[[i]])
jaccard_res(res_og_old_main2$row_clusters[[i]], res_og_old_main2$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_old_main_bis$row_clusters[[i]], res_old_main_bis$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(results$row_clusters[[i]], results$col_clusters[[i]], true_rows[[i]], true_cols[[i]])



jaccard_res(res_old_main$row_clusters[[i]], res_old_main$col_clusters[[i]], true_rows[[i]], true_cols[[i]])

i<- 1
jaccard_res(res_og_dg$row_clusters[[i]], res_og_dg$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_og_dg$Foutput[[i]]>1/200, res_og_dg$Goutput[[i]]>1/150, true_rows[[i]], true_cols[[i]])


jaccard_res(res_og40$row_clusters[[i]], res_og40$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_og50$row_clusters[[i]], res_og50$col_clusters[[i]], true_rows[[i]], true_cols[[i]])
jaccard_res(res_og40$Foutput[[i]]>1/200, res_og40$Goutput[[i]]>1/150, true_rows[[i]], true_cols[[i]])
jaccard_res(res_0$Foutput[[i]]>1/200, res_0$Goutput[[i]]>1/150, true_rows[[i]], true_cols[[i]])
jaccard_res(res_0$row_clusters[[i]], res_0$col_clusters[[i]], true_rows[[i]], true_cols[[i]])

res_0 <- restMultiNMTF_run(data, phi=0*phi_mat)

sils <- sil_score(data[[3]], res_og2$Foutput[[3]]>1/200,( res_og2$Goutput[[3]]>1/150)[,c(1,2,3,5,4)])
source("Functions/bisilhouette.r")
jaccard_res(res_og2$row_clusters[[2]], res_og2$col_clusters[[2]], true_rows[[2]], true_cols[[2]])
jaccard_res(res_og2$row_clusters[[3]], res_og2$col_clusters[[3]], true_rows[[3]], true_cols[[3]])
jaccard_res(res_og2$row_clusters[[1]], res_og2$col_clusters[[1]], true_rows[[1]], true_cols[[1]])

jaccard_res(res_og2$Foutput[[1]]>1/200, res_og2$Goutput[[1]]>1/250, true_rows[[1]], true_cols[[1]])
score <- sil_score(data[[1]], res_og2$row_clusters[[1]], res_og2$col_clusters[[1]])
score2 <- sil_score(data[[2]], res_og2$row_clusters[[2]], res_og2$col_clusters[[2]])
score3 <- sil_score(data[[3]], res_og2$row_clusters[[3]], res_og2$col_clusters[[3]])

j<- 3
i <- 5
max((data[[j]])[(1:200)[true_rows[[j]][,i]==1],(1:150)[true_cols[[j]][,i]==1]])

# generate synthetic data via our approach. 
<<<<<<< HEAD
source("SimStudy/RunSim/resNMTF_funcs.r")
source("resNMTF_funcs.r")

source("resNMTF_funcs.r")

=======
source("main.r")
>>>>>>> f653a3377b6b395d2c9bc5f87ccbfca8696c0bc3
source('SimStudy/RunSim/Functions/data_generation.r')
source('data_generation_og.r')
# generate a test dataset
n_views <- 5
row_cl_dims <- rep(200, n_views)
#100, 50,250 features respectively
col_cl_dims <- c(100, 50, 250)
<<<<<<< HEAD
col_cl_dims <- rep(200,3)
save_data(row_cl_dims, col_cl_dims, 5, 'test_data', 5, col_same_shuffle=FALSE)
data <- import_matrix("test_data/data.xlsx")
true_rows <- import_matrix("test_data/true_rows.xlsx")
true_cols <- import_matrix("test_data/true_cols.xlsx")
=======
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
>>>>>>> f653a3377b6b395d2c9bc5f87ccbfca8696c0bc3

max(data[[1]])
#perform gfa
source('SimStudy/RunSim/OtherMethods/gfa_funcs.r')
source("SimStudy/RunSim/Functions/extra_funcs.r")
test <- gfa_apply(data, 10)
test_high <- gfa_apply(data, 10, 0.9)
# check results
true_rows <- import_matrix("test_data/true_rows.xlsx")
true_cols <- import_matrix("test_data/true_cols.xlsx")
jaccard_res(test_high$row_clusters[[1]], test_high$col_clusters[[1]], true_rows[[1]], true_cols[[1]])
jaccard_res(test_high$row_clusters[[2]], test_high$col_clusters[[2]], true_rows[[2]], true_cols[[2]])
jaccard_res(test$row_clusters[[1]], test$col_clusters[[1]], true_rows[[1]], true_cols[[1]])
jaccard_res(test$row_clusters[[2]], test$col_clusters[[2]], true_rows[[2]], true_cols[[2]])
col_cl_dims <- list(rep(150, 2), rep(150, 3),
                rep(150, 4), rep(150, 5))


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




