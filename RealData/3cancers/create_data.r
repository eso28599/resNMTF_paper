# combine exp datasets
datasets <- c(
  "aml", "breast", "colon", "gbm", "kidney", "liver",
  "lung", "melanoma", "ovarian", "sarcoma"
)

genes_replace <- function(x) {
  x <- sub("?", "X.", x, fixed = TRUE)
  x <- sub("|", ".", x, fixed = TRUE)
  x <- gsub("-", ".", x, fixed = TRUE)
  return(x)
}

# aml
aml_exp <- read.table("aml/exp")
aml_genes <- row.names(aml_exp)
aml_methy <- read.table("aml/methy")
aml_m_features <- row.names(aml_methy)
aml_mirna <- read.table("aml/mirna")
aml_m_features_mirna <- row.names(aml_mirna)
aml_cols <- intersect(colnames(aml_mirna), colnames(aml_exp))
aml_exp <- aml_exp[, aml_cols]
aml_mirna <- aml_mirna[, aml_cols]

# breast
breast_exp <- read.table("breast/exp")
breast_genes_og <- row.names(breast_exp)
breast_genes <- row.names(breast_exp)
breast_genes <- genes_replace(breast_genes)
rownames(breast_exp) <- breast_genes
breast_methy <- read.table("breast/methy")
breast_m_features <- row.names(breast_methy)
breast_mirna <- read.table("breast/mirna")
breast_m_features_mirna <- row.names(breast_mirna)
breast_cols <- intersect(colnames(breast_mirna), colnames(breast_exp))
breast_exp <- breast_exp[, breast_cols]
breast_mirna <- breast_mirna[, breast_cols]
rownames(breast_mirna) <- genes_replace(row.names(breast_mirna))

# colon
colon_exp <- read.table("colon/exp")
colon_genes <- genes_replace(row.names(colon_exp))
colon_mirna <- read.table("colon/mirna")
rownames(colon_mirna) <- genes_replace(row.names(colon_mirna))
colon_m_features_mirna <- row.names(colon_mirna)
colon_cols <- intersect(colnames(colon_mirna), colnames(colon_exp))
colon_exp <- colon_exp[, colon_cols]
colon_mirna <- colon_mirna[, colon_cols]


# gene x individual
full_exp <- cbind(aml_exp, breast_exp, colon_exp)
zero_genes <- rowSums(full_exp)
selected_genes <- aml_genes[zero_genes > 0]
full_exp <- full_exp[selected_genes, ]

# mirna
selected_mirna <- intersect(aml_m_features_mirna, colon_m_features_mirna)
selected_mirna <- intersect(selected_mirna, breast_m_features_mirna)
zero_mirna <- rowSums(full_mirna)
selected_mirna <- selected_mirna[zero_mirna > 0]
full_mirna <- cbind(
  aml_mirna[selected_mirna, ],
  breast_mirna[selected_mirna, ],
  colon_mirna[selected_mirna, ]
)
openxlsx::write.xlsx(
  list(
    full_exp = full_exp,
    full_mirna = full_mirna
  ),
  file = "tgca_3cancers.xlsx",
  colNames = TRUE,
  rowNames = TRUE
)

true_labels <- matrix(0, nrow = ncol(full_exp), ncol = 3)
true_labels[1:ncol(aml_exp), 1] <- 1
true_labels[(ncol(aml_exp) + 1):(ncol(aml_exp) + ncol(breast_exp)), 2] <- 1
true_labels[(ncol(aml_exp) + ncol(breast_exp) + 1):ncol(full_exp), 3] <- 1
write.csv(
  true_labels,
  file = "tgca_3cancers_labels.csv",
  row.names = FALSE
)




# length(selected_mirna)

# colnames(breast_mirna)
# all(aml_genes == colon_genes)

# dim(full_exp)
# # check for zero genes


# full_mat <- unname(as.matrix(full_exp))

# gene_sd <- apply(full_exp, 1, sd)
# gene_mean <- apply(full_exp, 1, mean)
# sum(gene_sd / gene_mean < 1)

# library(biclust)

# test <- biclust(full_mat, method = BCPlaid)

# library(resnmtf)

# nmtf_results <- apply_resnmtf(list(full_mat), k_vec = c(3), stability = FALSE)


# # gbm
# gbm_exp <- read.table("gbm/exp")
# length(intersect(aml_genes, gmb_genes))
# length(gmb_genes) # 12042
# gmb_genes <- genes_replace(row.names(gbm_exp))
# length(intersect(aml_genes, gmb_genes)) == length(gmb_genes)

# # kidney
# kidney_exp <- read.table("kidney/exp")
