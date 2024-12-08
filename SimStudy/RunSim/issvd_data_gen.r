args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
batch_folder = as.character(args[2])
case = as.character(args[3])
source(paste0(path_to_sim_folder,
             "/sim_parameters.r")) #parameters for simulation
library(openxlsx)
library(MASS)
#iSSVD scenario 1
#X1 100 x 1000
#X2 100 x 1000
#4 biclusters
#each bicluster is 10x100

#U is 100x100
#we randomly select 10 rows in each column of the first four columns in matrix 𝑈
# and assign data values generated from a uniform distribution 𝑈(0.5,1)
#, while ensuring there is no overlapping samples.
true_rows <- matrix(0, 100, 4)
true_cols <-list(matrix(0, 1000, 4), matrix(0, 1000, 4))
U <- matrix(0, 100, 100)
full_list <- 1:100
biclustered <- c(150)
for(i in 1:4){
    non_zeros <- sample(full_list[-biclustered], 10)
    U[non_zeros, i] <- runif(10, 0.5, 1)
    biclustered <- c(biclustered, non_zeros)
    true_rows[non_zeros, i] <- 1
}
U[, 5:100] <- mvrnorm(100, rep(0, 96), diag(rep(1, 96)))
true_rows <- list(true_rows, true_rows)

#1000 x 100
V1 <- matrix(0, 1000, 100)
full_list <- 1:1000
biclustered <- c(1500)
for(i in 1:4){
    non_zeros <- sample(full_list[-biclustered], 100)
    V1[non_zeros, i] <- runif(100, 0.5, 1)
    biclustered <- c(biclustered, non_zeros)
    true_cols[[1]][non_zeros, i] <- 1
}
V1[, 5:100] <- mvrnorm(1000, rep(0, 96), diag(rep(1, 96)))

#1000 x 100
V2 <- matrix(0, 1000, 100)
full_list <- 1:1000
biclustered <- c(1500)
for(i in 1:4){
    non_zeros <- sample(full_list[-biclustered], 100)
    V2[non_zeros, i] <- runif(100, 0.5, 1)
    biclustered <- c(biclustered, non_zeros)
    true_cols[[2]][non_zeros, i] <- 1
}
V2[, 5:100] <- mvrnorm(1000, rep(0, 96), diag(rep(1, 96)))

#singular values
eps <- 0.3
S <- diag(c(27, 20, 18, 10, rep(eps, 96)))

X1 <- U %*% S %*% t(V1)
if(case=="two"){
    X2 <- 5*(U %*% S %*% t(V2))
    #save data
    for(i in 1:length(sigma_vec)){
        #noise
        sigma2 <- sigma_vec[i]^2
        noise1 <- mvrnorm(100, rep(0, 1000), sigma2*diag(1, 1000))
        noise2 <- mvrnorm(100, rep(0, 1000), sigma2*diag(1, 1000))
        data_mod <- list("v1"=(X1 + noise1), "v2"=(X2 + noise2))
        file_path <- paste0(path_to_save, method_vec[i])
        openxlsx::write.xlsx(data_mod, file = paste0(file_path, "/data.xlsx"))
        openxlsx::write.xlsx(true_rows, file = paste0(file_path, "/true_rows.xlsx")) 
        openxlsx::write.xlsx(true_cols, file = paste0(file_path, "/true_cols.xlsx"))
    }
}else{
    sigma2 <- 0.2^2
    noise1 <- mvrnorm(100, rep(0, 1000), sigma2*diag(1, 1000))
    noise2 <- mvrnorm(100, rep(0, 1000), sigma2*diag(1, 1000))
    X2 <- (U %*% S %*% t(V2))
    for(i in 1:length(scalar_vec)){
        data_mod <- list("v1"=(X1 + noise1), "v2"=(scalar_vec[i]*X2 + noise2))
        file_path <- paste0(path_to_save, method_vec[i])
        openxlsx::write.xlsx(data_mod, file = paste0(file_path, "/data.xlsx"))
        openxlsx::write.xlsx(true_rows, file = paste0(file_path, "/true_rows.xlsx")) 
        openxlsx::write.xlsx(true_cols, file = paste0(file_path, "/true_cols.xlsx"))
    }
}
