# combine results
filepath <- "RealDate/3sources"
resnmtf_res <- read.csv(paste0(filepath,
 "/stability_study.csv"))[,c(1:16,21:23)]
gfa_res <- read.csv(paste0(filepath, "/gfa_method_results.csv"))[1:5, ]
issvd_res <- read.csv(paste0(filepath, "/python_method_results.csv"))

colnames <- c(
  "Method",
  paste0("F score (V", 1:n_views, ")"), "F score",
  paste0("Relevance (V", 1:n_views, ")"), "Relevance",
  paste0("Recovery (V", 1:n_views, ")"), "Recovery",
  paste0("BiS (V", 1:n_views, ")"), "BiS"
)

gfa_result <- c("GFA", as.matrix(gfa_res[
  which.max(gfa_res$F.score),
  4:ncol(gfa_res)
]))
issvd_result <- c("iSSVD", as.matrix(issvd_res[
  which.max(issvd_res$F.score),
  4:ncol(issvd_res)
]))
resnmtf_results <- as.matrix(resnmtf_res[, 2:ncol(resnmtf_res)])

all_results <- rbind(
  resnmtf_results, conc_result,
  gfa_result, issvd_result
)
colnames(all_results) <- colnames
all_results <- as.data.frame(all_results)
all_results[, 2:ncol(all_results)] <- apply(
  all_results[, 2:ncol(all_results)],
  2, as.numeric
)
write.csv(all_results, paste0(filepath, "/all_results.csv"), row.names = FALSE)
res_text <- kbl(all_results,
  booktabs = TRUE, "latex", escape = FALSE,
  digits = 4, row.names = FALSE
)
df <- kable_styling(res_text)
# save
sink(paste0(filepath, "/all_results.txt"))
print(df)
sink()

colnames_sub <- c("Method", "F score", "Relevance", "Recovery", "BiS")
all_results_main <- all_results[, colnames_sub]
write.csv(all_results_main, paste0(filepath, "/all_results_main.csv"),
  row.names = FALSE
)
res_text <- kbl(all_results_main,
  booktabs = TRUE, "latex", escape = FALSE,
  digits = 4, row.names = FALSE
)
df <- kable_styling(res_text)
# save
sink(paste0(filepath, "/all_results_main.txt"))
print(df)
sink()

colnames_sub <- c("Method", paste0("F score (V", 1:n_views, ")"), "F score")
all_results_main <- all_results[, colnames_sub]
write.csv(all_results_main, paste0(filepath, "/all_results_fscores.csv"),
  row.names = FALSE
)
res_text <- kbl(all_results_main,
  booktabs = T, "latex", escape = FALSE,
  digits = 4, row.names = FALSE
)
df <- kable_styling(res_text)
# save
sink(paste0(filepath, "/all_results_fscores.txt"))
print(df)
sink()