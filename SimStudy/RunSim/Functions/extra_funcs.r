# additional functions and packages
library(curl)
library(readxl)
library(rio)
library(dplyr)
library(kableExtra)
library(latex2exp)
library(ggpubr)
library(philentropy)

# import as list of matrices instead of as a list
import_matrix <- function(filename, col_n = TRUE) {
    return(lapply(
        import_list(filename, col_names = col_n),
        function(x) as.matrix(x)
    ))
}
save_results <- function(results, file_path, error = T) {
    row_filename <- paste0(file_path, "/row_clusts.xlsx")
    col_filename <- paste0(file_path, "/col_clusts.xlsx")
    error_filename <- paste0(file_path, "/errors.csv")
    openxlsx::write.xlsx(results$row_clusters, file = row_filename)
    openxlsx::write.xlsx(results$col_clusters, file = col_filename)
    if (error) {
        suppressMessages(write.csv(results$Error, file = error_filename))
    }
}

save_results_real <- function(results, file_path, index, error = T) {
    row_filename <- paste0(file_path, "/row_clusts_", index, ".xlsx")
    col_filename <- paste0(file_path, "/col_clusts_", index, ".xlsx")
    openxlsx::write.xlsx(results$row_clusters, file = row_filename)
    openxlsx::write.xlsx(results$col_clusters, file = col_filename)
}
