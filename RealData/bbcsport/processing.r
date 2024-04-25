#bbc processing
library(R.matlab)
bbc <- readMat("investigate_applications/BBC/bbcsport_2view.mat")
bbc_data <- list(bbc$data[[1]][[1]], bbc$data[[2]][[1]])
zeros_cols_1 <- (1:3183)[(colSums(bbc_data[[1]])!=0)]
zeros_cols_2 <- (1:3203)[(colSums(bbc_data[[2]])!=0)]
bbc_data[[1]] <- bbc_data[[1]][,zeros_cols_1]
bbc_data[[2]] <- bbc_data[[2]][,zeros_cols_2]
bbc_d2 <- list(t(bbc_data[[1]]), t(bbc_data[[2]]))
bbc_d2 <- lapply(bbc_d2, as.matrix)
bbc_rows <- matrix(0, 544, 5)
for(i in 1:5){
  bbc_rows[,i] <- ifelse(bbc$truth==i, 1, 0)
}
write.csv(bbc_rows, "investigate_applications/BBC/bbc_rows_truth.csv")
openxlsx::write.xlsx(bbc_d2, file = paste0("investigate_applications/BBC/", "bbc_data_processed.xlsx"))
