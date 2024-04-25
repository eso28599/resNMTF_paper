#convergence results
source("sim_common_files/resNMTF_funcs.r")
source("sim_common_files/data_generation.r")
source("sim_common_files/extra_funcs.r")
row_start <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end <- list(c(75, 150, 200), c(75, 150, 200), c(75, 150, 200))
col_start <- list(c(1, 31, 61), c(1, 11, 31), c(1, 101, 151))
col_end <- list(c(30, 60, 100), c(10, 30, 50), c(100, 150, 250))

phi <- matrix(0, nrow = 3, ncol = 3)
phi[1, c(2, 3)] <- 200
phi[2, c(3)] <- 200

save_data(rep(200,3), c(100, 50, 250),
                 row_start, row_end,
                  col_start, col_end,
                "results", 2,  TRUE, FALSE)

data <- import_matrix("results/data.xlsx")
nmtf_results <- restMultiNMTF_run(Xinput = data, phi = phi)
evaluate_simulation_comp(nmtf_results$row_clusters, nmtf_results$col_clusters,
                                                import_matrix("results/true_rows.xlsx"), import_matrix("results/true_cols.xlsx"), data)
row_start_5 <- list(c(1, 51, 76, 111, 151), c(1, 51, 76, 111, 151),
                     c(1, 51, 76, 111, 151))
row_end_5 <- list(c(50, 75, 110, 150, 200), c(50, 75, 110, 150, 200),
                 c(50, 75, 110, 150, 200))
col_start_5 <- list(c(1, 21, 31, 46, 61), c(1, 6, 11, 21, 31),
                         c(1, 51, 101, 126, 151))
col_end_5 <- list(c(20, 30, 45, 60, 100), c(5, 10, 20, 30, 50),
                     c(50, 100, 125, 150, 250))

save_data(rep(200,3), c(100, 50, 250),
                 5,
                "results/5", 2,  TRUE, FALSE)
save_data(rep(200,3), c(100, 50, 250),
                 row_start_5, row_end_5,
                  col_start_5, col_end_5,
                "results/5", 2,  TRUE, FALSE)

data5 <- import_matrix("results/5/data.xlsx")
nmtf_results5 <- restMultiNMTF_run(Xinput = data5, phi = 0.5*phi)

source("sim_common_files/evaluation_funcs.r")
evaluate_simulation_comp(nmtf_results5$row_clusters, nmtf_results5$col_clusters,
                                                import_matrix("results/5/true_rows.xlsx"), import_matrix("results/5/true_cols.xlsx"), data5)
row_start_6 <- list(c(1, 51, 76, 111, 151, 176), c(1, 51, 76, 111, 151, 176),
                     c(1, 51, 76, 111, 151, 176))
row_end_6 <- list(c(50, 75, 110, 150, 175, 200), c(50, 75, 110, 150, 175, 200),
                 c(50, 75, 110, 150, 175, 200))
col_start_6 <- list(c(1, 21, 31, 46, 61, 81), c(1, 6, 11, 21, 31, 41),
                         c(1, 51, 101, 126, 151, 201))
col_end_6 <- list(c(20, 30, 45, 60, 80, 100), c(5, 10, 20, 30, 40, 50),
                     c(50, 100, 125, 150, 200, 250))    

save_data(rep(200,3), c(100, 50, 250),
                 row_start_6, row_end_6,
                  col_start_6, col_end_6,
                "results/6", 2,  TRUE, FALSE)

data6 <- import_matrix("results/6/data.xlsx")

source("sim_common_files/evaluation_funcs.r")
evaluate_simulation_comp(nmtf_results6$row_clusters, nmtf_results6$col_clusters,
                                                import_matrix("results/6/true_rows.xlsx"), import_matrix("results/6/true_cols.xlsx"), data6)  
nmtf_res <- vector("list", length=4)
i<-1
for(phi_val in c(100, 300, 500, 700)){
  nmtf_res[[i]] <- restMultiNMTF_run(Xinput = data6, phi = phi_val/200*phi, nIter = 200)
  i<- i+1
  evaluate_simulation_comp(nmtf_res[[i]]$row_clusters, nmtf_res[[i]]$col_clusters,
                                                import_matrix("results/6/true_rows.xlsx"), import_matrix("results/6/true_cols.xlsx"), data6)  
}
nmtf_results6_0 <- restMultiNMTF_run(Xinput = data6, phi = phi, nIter = 200)
source("sim_common_files/evaluation_funcs.r")
evaluate_simulation_comp(nmtf_results6_0$row_clusters, nmtf_results6_0$col_clusters,
                                                import_matrix("results/6/true_rows.xlsx"), import_matrix("results/6/true_cols.xlsx"), data6)                                     

row_start_2 <- list(c(1, 76), c(1, 76), c(1, 76))
row_end_2 <- list(c(75, 200), c(75, 200), c(75, 200))
col_start_2 <- list(c(1, 31), c(1, 11), c(1, 101))
col_end_2 <- list(c(30, 100), c(10, 50), c(100, 250))
save_data(rep(200,3), c(100, 50, 250),
                 row_start_2, row_end_2,
                  col_start_2, col_end_2,
                "results/2b", 2,  TRUE, FALSE)

data2 <- import_matrix("results/2/data.xlsx")
nmtf_results2 <- restMultiNMTF_run(Xinput = data2, phi = phi, k_min=2, nIter=200)
evaluate_simulation_comp(nmtf_results2$row_clusters, nmtf_results2$col_clusters,
                                                import_matrix("results/2/true_rows.xlsx"), import_matrix("results/2/true_cols.xlsx"), data2)   

data2b <- import_matrix("results/2b/data.xlsx")
nmtf_results2b <- restMultiNMTF_run(Xinput = data2b, phi = phi, k_min=3, nIter=200)
evaluate_simulation_comp(nmtf_results2b$row_clusters, nmtf_results2b$col_clusters,
                                                import_matrix("results/2b/true_rows.xlsx"), import_matrix("results/2b/true_cols.xlsx"), data2b)         

row_start_4 <- list(c(1, 76, 111, 151), c(1, 76, 111, 151), c(1, 76, 111, 151))
row_end_4 <- list(c(75, 110, 150, 200), c(75, 110, 150, 200),
                 c(75, 110, 150, 200))
col_start_4 <- list(c(1, 31, 46, 61), c(1, 11, 21, 31), c(1, 101, 126, 151))
col_end_4 <- list(c(30, 45, 60, 100), c(10, 20, 30, 50), c(100, 125, 150, 250))
save_data(rep(200,3), c(100, 50, 250),
                 row_start_4, row_end_4,
                  col_start_4, col_end_4,
                "results/4", 4,  TRUE, FALSE)

data4 <- import_matrix("results/4/data.xlsx")
nmtf_results4 <- restMultiNMTF_run(Xinput = data4, phi = phi)
evaluate_simulation_comp(nmtf_results4$row_clusters, nmtf_results4$col_clusters,
                                                import_matrix("results/4/true_rows.xlsx"), import_matrix("results/4/true_cols.xlsx"), data4)      


library(ggplot2)
n_its<- 1000
errors <- c(nmtf_results2$All_Error[1:n_its], nmtf_results$All_Error[1:n_its],  nmtf_results4$All_Error[1:n_its], 
nmtf_results5$All_Error[1:n_its], nmtf_results6$All_Error[1:n_its])
bicls <- factor(c(rep(2, n_its), rep(3, n_its), rep(4, n_its), 
        rep(5, n_its), rep(6, n_its)))
its <- rep(1:n_its, 5)
df <-  data.frame("Error"=errors, "Biclusters" = bicls, "Iteration"=its)
write.csv(df, "results/error_df.csv")
df <- read.csv("results/error_df.csv" )
df$Biclusters <- factor(df$Biclusters)
df2 <- subset(df, Iteration<=250)
p <- ggplot(data=df2, aes(x=Iteration, y=Error, group=Biclusters)) +
    geom_line(aes(color=Biclusters), linewidth=0.5)+ggtitle("Error across iterations")

p <- ggplot(data=df2, aes(x=Iteration, y=Error, group=Biclusters)) +
    geom_line(aes(color=Biclusters), linewidth=0.5)
ggsave("results/error_plot.pdf", plot = p, compress = FALSE, device="pdf", width=5, height=5)



new_sil_score <- function(Xinput, Foutput, Goutput, row_clustering, col_clustering, index){
    #simultaneously calculates silhouette score for a 
    #clustering as well matching clusters correctly.
    n_views <- length(Xinput)
    n_clusts <- dim(Foutput[[1]])[2]
    relations2 <- matrix(0, nrow = n_views, n_clusts)
    relations3 <- matrix(0, nrow = n_views, n_clusts)
    sil_score2 <- matrix(0, nrow = n_views, n_clusts)
    sil_score3 <- matrix(0, nrow = n_views, n_clusts)
    if (index == 2) {
        clust_one <- col_clustering
        clust_two <- row_clustering
    }else{
        clust_one <- row_clustering
        clust_two <- col_clustering
    }
    for (i in 1:n_views){
        s_mat <- matrix(0, nrow = n_clusts, n_clusts)
        a_mat <- matrix(0, nrow = n_clusts, n_clusts)
        b_mat <- matrix(0, nrow = n_clusts, n_clusts)
        for (k in 1:n_clusts){
          #select data from specific column clustering
          if (sum(clust_one[[i]][, k] == 1) == 0){
            #print("no cols")
            s_mat[k, ] <- 0
          }else{
            if (index == 2){
              new_data <- Xinput[[i]][, (clust_one[[i]][, k] == 1)]
            }else{
              new_data <- t(Xinput[[i]][(clust_one[[i]][, k] == 1), ])
            }
            spear_dists <- as.matrix(dist(new_data, upper = FALSE, diag = TRUE))/ (dim(new_data)[2])
            for (j in 1:n_clusts){
              indices <- clust_two[[i]][, j] == 1
              if (sum(indices) == 0) {
                #("none")
                s_mat[k, j] <- 0
              }else{
                #a_vals <- rowMeans(spear_dists[indices, indices])
                a_vals <- apply(spear_dists[indices, indices], 1, function(x) sum(x)/(length(x)-1))
                #other clusts
                other <- (1:n_clusts)[-j]
                b_vec <- c()
                b_vals <- vector("list", length = (n_clusts - 1))
                t <- 1
                for(l in other){
                    oth_ind <- clust_two[[i]][, l] == 1
                    b_val <- rowMeans(spear_dists[indices, oth_ind])
                    b_vec <- c(b_vec, mean(b_val))
                    b_vals[[t]] <- b_val
                    t <- t + 1
                }
                closest <- which.min(b_vec)
                #print(b_vec)
                b_vals <- b_vals[[closest]]
                s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
                s_mat[k, j] <- mean(s_con)
                #a_mat[k, j] <- mean(a_vals)
                #b_mat[k,j] <- mean(b_vals)
               }
            }
          }
        }
        #print(s_mat)
        if(n_clusts == 2){
          relations[i, ] <- apply(s_mat, 1, which.min)
          sil_score[i, ] <- apply(s_mat, 1, min)
        }else{
          relations[i, ] <- apply(s_mat, 1, which.max)
          sil_score[i, ] <- apply(s_mat, 1, max)
        }
        
      
      #find the contribution to the s_score of the biclusters from this view
      
    }
  acc_vals <- sil_score[sil_score!=0]
  sil <- ifelse(n_clusts == 1, mean(acc_vals),
                 mean(acc_vals) - 2 * sd(acc_vals))
  #sil <- ifelse(n_clusts == 1, mean(sil_score[sil_score!=0]),
                 #mean(sil_score) - 2 * sd(sil_score))
  #return relationships and mean of sil_score across views
  return(list("sil" = sil, "scores" = sil_score, "relations" = relations))
}
