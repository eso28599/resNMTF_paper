#single cell omics from theo's thesis
load("investigate_applications/single_cell/data_A549.rda")
source("sim_common_files/resNMTF_files/resNMTF_funcs.r")
source("sim_common_files/resNMTF_files/evaluation_funcs.r")
v1 <- data_A549[[1]][[1]]
v2 <- data_A549[[1]][[2]]
labs <- data_A549[[2]]

#remove those with low counts
low_count1 <- (rowSums(v1!=0)/(dim(v1)[1]))<= 0.05
v1 <- v1[!low_count1, ]
lab_mat <- cbind(ifelse(labs==0,1,0), 
                    ifelse(labs==1,1,0),
                    ifelse(labs==3,1,0))
sd <- apply(v2, 1, sd)
v2_small <- v2[sd>0.4,]
#save

openxlsx::write.xlsx(list(t(as.matrix(v1)), t(as.matrix(v2_small))), 
            file = "investigate_applications/single_cell/data_processed.xlsx")
write.csv(lab_mat,"investigate_applications/single_cell/true_labs.csv")

labs <- read.csv("investigate_applications/single_cell/true_labs.csv")[,2:4]
#remove those individuals with low counts
#plot(colSums(v2))
#min(colSums(v2))
phi_mat <- matrix(0,2,2)
phi_mat[1,2] <- 1
#reduce <- prcomp(t(v2_small))
v1 <- matrix(v1)
test1 <- restMultiNMTF_run(list(t(v1), t(v2_small)), k_min=3, k_max=4, phi=700*phi_mat, stability=FALSE)
test1000_diff_data <- restMultiNMTF_run(list(t(v1), t(v2_small)), k_min=4, k_max=4, phi=1000*phi_mat, stability=FALSE)

test1000 <- restMultiNMTF_run(list(t(v1), t(v2_small)), k_min=4, k_max=4, phi=1000*phi_mat, stability=FALSE, distance="manhattan")
test1_0 <- restMultiNMTF_run(list(t(v1), t(v2_small)), k_min=3, k_max=4, phi=0*phi_mat, stability=FALSE)
test1 <- restMultiNMTF_run(list(t(v2_small)), k_min=3, k_max=4, stab_thres=0.05)

jaccard_row(test1000$Foutput[[2]]>1/2641, lab_mat)
sil_score(t(as.matrix(v1)), test1000$row_clusters[[1]], test1000$col_clusters[[1]], method="euclidean")$sil

jaccard_row(test1000$row_clusters[[1]], lab_mat)
test1_t <- restMultiNMTF_run(list(v1, (v2_small)), k_min=3, k_max=4, psi=0*phi_mat, stability=FALSE)
sil1 <- sil_score(as.matrix(v1), test1$col_clusters[[1]], test1$row_clusters[[1]], method="euclidean")
sil1 <- sil_score(as.matrix(v1), test1$col_clusters[[1]], test1$row_clusters[[1]], method="manhattan")
sil_score(t(as.matrix(v1)),  test1$row_clusters[[1]], test1$col_clusters[[1]], method="manhattan")$sil
sil1t <- sil_score(as.matrix(v1), test1_t$row_clusters[[1]], test1_t$col_clusters[[1]], method="manhattan")
library(ggfortify)
dim(v1)
#SeuratData
install.packages("SeuratData")
library(SeuratData)


