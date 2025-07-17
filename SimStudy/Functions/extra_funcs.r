library(openxlsx)
library(rio)

# import as list of matrices instead of as a list
import_matrix <- function(filename, col_n = TRUE) {
  return(lapply(
    import_list(filename, col_names = col_n), function(x) as.matrix(x)
  ))
}

import_matrix <- function(csv_dir, base_name, row_names = NULL) {
  # from chatgpt
  # List all CSV files that match the base name pattern
  csv_files <- list.files(path = csv_dir, pattern = paste0("^", base_name, "_.*\\.csv$"), full.names = TRUE)
  
  # Initialize an empty list
  matrix_list <- list()
  
  # Loop over each CSV file and read it as a matrix
  for (file in csv_files) {
    # Extract sheet name from filename
    sheet_name <- sub(paste0("^", base_name, "_"), "", basename(file))
    sheet_name <- sub("\\.csv$", "", sheet_name)
    
    # Read the CSV and convert to matrix
    df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE, row.names = row_names)
    matrix_list[[sheet_name]] <- as.matrix(df)
  }
  
  return(matrix_list)
}


export_matrix_list <- function(mat_list, base_name, output_dir = ".") {
  for (sheet_name in seq_along(mat_list)) {
    filename <- file.path(output_dir, paste0(base_name, "_", sheet_name, ".csv"))
    # Write matrix to CSV
    write.csv(mat_list[[sheet_name]], file = filename, row.names = FALSE)
    message("Saved matrix: ", filename)
  }
}


save_results <- function(results, file_path, error = TRUE) {
  row_filename <- paste0(file_path, "/row_clusts.xlsx")
  col_filename <- paste0(file_path, "/col_clusts.xlsx")
  error_filename <- paste0(file_path, "/errors.csv")
  openxlsx::write.xlsx(results$row_clusters, file = row_filename)
  openxlsx::write.xlsx(results$col_clusters, file = col_filename)
  if (error) {
    suppressMessages(write.csv(results$Error, file = error_filename))
  }
}

save_results_real <- function(results, file_path, index, error = TRUE) {
  row_filename <- paste0(file_path, "/row_clusts_", index, ".xlsx")
  col_filename <- paste0(file_path, "/col_clusts_", index, ".xlsx")
  openxlsx::write.xlsx(results$row_clusters, file = row_filename)
  openxlsx::write.xlsx(results$col_clusters, file = col_filename)
}