# read txt file
data <- read.table("proteomics.txt")
# remove the first column
protein_names <- data[,1]
identity <- data[1,]
protein_data <- apply(data[2:nrow(data),2:ncol(data)],2,as.numeric)

heatmap(protein_data[1:50,], Rowv=NA, Colv=NA, col=cm.colors(256), scale="column")
plot(apply(protein_data,1, mean))

mean <- apply(protein_data,1, mean)

# apply resnmtf
source("main.r")
protein_data <- protein_data[,1:436]
res_t <- restMultiNMTF_run(list(t(protein_data)),repeats = 5, no_clusts = FALSE, stability = FALSE)
res_og <- restMultiNMTF_run(list(protein_data),repeats = 5, no_clusts = FALSE, stability = FALSE)


partic <- read.csv("aals_participants.csv")
identity[1:5]
# extract everything after a _ 
identity <- gsub(".*_", "", identity)[1:5]
sub("CASE_", "", identity)[1:5]
(gsub(".*_", "", (identity)))[1:5]
#reverese identity
identity <- rev(identity)
rev(identity[1])[1]
id <- lapply(identity[2:437], function(x) strsplit(x, "_")[[1]][2])
id <- lapply(identity[2:437], function(x) x[1:4])
id <- lapply(id, function(x) strsplit(x, "-")[[1]][1])
gsub("_*.","",identity[2:437][1:5])
#extract first four characters of a string
id <- unlist(lapply(identity[2:437], function(x) substr(x, 1, 4)))
id <- ifelse(id == "CASE", 1, 0)
source("SimStudy/Functions/evaluation_funcs.r")
jaccard_row(res_og$col_clusters[[1]], cbind(id,1-id), print=TRUE)
jaccard_row(res_t$row_clusters[[1]], cbind(id,1-id), print=TRUE)

sum((rowSums(res$col_clusters[[1]])>1)==(1-id))
sum(res_og$col_clusters[[1]][,2]==res_og$col_clusters[[1]][,3])
sum((res$col_clusters[[1]][,2]==1 )& (res$col_clusters[[1]][,3]==1))
sum((res_og$col_clusters[[1]][,2]==1 )& (res_og$col_clusters[[1]][,3]==1))


plot(res$Goutput[[1]][,1])


source('SimStudy/Functions/data_generation.r')
# generate a test dataset
n_views <- 3
row_cl_dims <- rep(200, n_views)
#100, 50,250 features respectively
col_cl_dims <- c(100, 50, 250)
save_data(row_cl_dims, col_cl_dims, 5, 'test_data', 5 ,col_same_shuffle=FALSE)


source('main.r')
source('SimStudy/Functions/extra_funcs.r')
data <- import_matrix("test_data/data.xlsx")
test_init <- init_mats_random(data, 5)
