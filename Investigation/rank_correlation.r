# calculate rank correlation between fscore and bis for all scenarios
# for figure investigating base scenario of 3 views 5 biclusters

# bicl
res <- read.csv("SimStudy/RunSim/Results/bicl/bicl_3v/all_results.csv")
f_score <- res[res$Measure == "F_score",]
bis <- res[res$Measure == "BiS",]
rank_f_score <- c()
rank_bis <- c()
top_rank <- c()
for(k in c(3,4,5,6)){
    f_score_k <- f_score[f_score$Biclusters==k,]
    bis_k <- bis[bis$Biclusters==k,]
    rank_f_score <- c(rank_f_score, rank(f_score_k$Mean))
    rank_bis <- c(rank_bis, rank(bis_k$Mean))
    top_rank <- c(top_rank, which.max(f_score_k$Mean)==which.max(bis_k$Mean))
}

#views
res <- read.csv("SimStudy/RunSim/Results/views/views_5b/all_results.csv")
f_score <- res[res$Measure == "F_score",]
bis <- res[res$Measure == "BiS",]
for(k in c(2,3,4,5)){
    f_score_k <- f_score[f_score$Views==k,]
    bis_k <- bis[bis$Views==k,]
    rank_f_score <- c(rank_f_score, rank(f_score_k$Mean))
    rank_bis <- c(rank_bis, rank(bis_k$Mean))
    top_rank <- c(top_rank, which.max(f_score_k$Mean)==which.max(bis_k$Mean))
}

#indiv
res <- read.csv("SimStudy/RunSim/Results/indiv/indiv_3v5b/all_results.csv")
f_score <- res[res$Measure == "F_score",]
bis <- res[res$Measure == "BiS",]
for(k in c(50,200,300,500,1000)){
    f_score_k <- f_score[f_score$Views==k,]
    bis_k <- bis[bis$Views==k,]
    rank_f_score <- c(rank_f_score, rank(f_score_k$Mean))
    rank_bis <- c(rank_bis, rank(bis_k$Mean))
    top_rank <- c(top_rank, which.max(f_score_k$Mean)==which.max(bis_k$Mean))
}

# noise
res <- read.csv("SimStudy/RunSim/Results/noise/noise_3v5b/all_results.csv")
f_score <- res[res$Measure == "F_score",]
bis <- res[res$Measure == "BiS",]
for(k in c(1:9, 10*(1:10))){
    f_score_k <- f_score[f_score$Noise==k,]
    bis_k <- bis[bis$Noise==k,]
    rank_f_score <- c(rank_f_score, rank(f_score_k$Mean))
    rank_bis <- c(rank_bis, rank(bis_k$Mean))
    top_rank <- c(top_rank, which.max(f_score_k$Mean)==which.max(bis_k$Mean))
}




# correlation plots
bicl_cor <- read.csv("SimStudy/RunSim/Results/bicl/bicl_3v/biclcorr_tab.txt")
bicl_f <- as.numeric(substr(bicl_cor[4, ], 11, 19))
bicl_rel <- as.numeric(substr(bicl_cor[5, ], 13, 21))

#views
views_cor <- read.csv("SimStudy/RunSim/Results/views/views_5b/viewscorr_tab.txt")
views_f <- as.numeric(substr(views_cor[4, ], 11, 19))
views_rel <- as.numeric(substr(views_cor[5, ], 13, 21))

#indiv
indiv_cor <- read.csv("SimStudy/RunSim/Results/indiv/indiv_3v5b/biclcorr_tab.txt")
indiv_f <- as.numeric(substr(indiv_cor[4, ], 11, 19))
indiv_rel <- as.numeric(substr(indiv_cor[5, ], 13, 21))

# noise
noise_cor <- read.csv("SimStudy/RunSim/Results/noise/noise_3v5b/noisecorr_tab.txt")
noise_f <- as.numeric(substr(noise_cor[4, ], 11, 19))
noise_rel <- as.numeric(substr(noise_cor[5, ], 13, 21))

(noise_f + bicl_f + views_f + indiv_f)/4
(noise_rel + bicl_rel + views_rel + indiv_rel)/4

# sc <- import_matrix("RealData/single_cell/data_processed.xlsx")
# heatmap(sc[[1]])
# cor(sc[[1]])
