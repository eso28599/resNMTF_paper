#3sources data processing
library(Matrix)

load_three <- function(outlet){
    name <- paste0("3sources/3sources_", outlet)
    view <- readMM(paste0(name, ".mtx"))
    docs <- read.csv(paste0(name, ".docs"), header=FALSE)
    terms <- read.csv(paste0(name, ".terms"), header=FALSE)  
    return(list("view"=view, "docs"=docs, "terms"=terms))
}

bbc <- load_three("bbc")
guardian <- load_three("guardian")
reuters <- load_three("reuters")

docs <- intersect(intersect(reuters$docs[,1], guardian$docs[,1]), bbc$docs[,1])
terms <- intersect(intersect(reuters$terms[,1], guardian$terms[,1]), bbc$terms[,1])

update_three <- function(data, doc_vec, term_vec){
    data$docs <- apply(data$docs, 1, function(x) ifelse(any(x==(doc_vec)), TRUE, FALSE))
    data$terms <- apply(data$terms, 1, function(x) ifelse(any(x==(term_vec)), TRUE, FALSE))
    data$view <- matrix(data$view, dim(data$view)[1], dim(data$view)[2])[data$terms, data$docs]
    return(data)
}

bbc <- update_three(bbc, docs, terms)
guardian <- update_three(guardian, docs, terms)
reuters <- update_three(reuters, docs, terms)

#find which terms any of the views have not used
zeros <-(1:length(terms))[(rowSums(bbc$view)==0)|(rowSums(reuters$view)==0)|(rowSums(guardian$view)==0)]
three_source <- list("b" = bbc$view[-zeros,], "g" = guardian$view[-zeros,], "r" = reuters$view[-zeros,])
#save data
openxlsx::write.xlsx(three_source, "3sources/3sources_all.xlsx")

#true labels
overlap_labels <- t(read.csv("3sources/3sources.overlap.clist", header=FALSE))
overlap_labels[1, 1:6] <- c(1, 5, 42, 1, 6, 2)

find_members <- function(labels, i){
        intersect(as.numeric(labels[,i][!is.na(labels[,i])]), docs)
}
c_vecs <- vector("list", length=6)
true_labs <- matrix(0, length(docs), 6) 
for(i in 1:6){
    c_vecs[[i]] <- find_members(overlap_labels, i)
    true_labs[sapply(docs,  function(x) ifelse(any(c_vecs[[i]]==(x)), TRUE, FALSE)),i] <- 1 
}

write.csv(true_labs, "3sources/true_labels.csv")
source("ResNMTF/ResNMTF/SimStudy/Functions/extra_funcs.r")
sources <- import_matrix("ResNMTF/ResNMTF/RealData/3sources/3sources_all_diff.xlsx")
s1 <- sum(sources[[1]]==0)/(dim(sources[[1]])[1]*dim(sources[[1]])[2])
s2 <- sum(sources[[2]]==0)/(dim(sources[[2]])[1]*dim(sources[[2]])[2])
s3 <- sum(sources[[3]]==0)/(dim(sources[[3]])[1]*dim(sources[[3]])[2])
bbc <- import_matrix("ResNMTF/ResNMTF/RealData/bbcsport/bbc_data_processed.xlsx")
b1 <- sum(bbc[[1]]==0)/(dim(bbc[[1]])[1]*dim(bbc[[1]])[2])
b2 <- sum(bbc[[2]]==0)/(dim(bbc[[2]])[1]*dim(bbc[[2]])[2])
a549 <- import_matrix("ResNMTF/ResNMTF/RealData/single_cell/data_processed.xlsx")
a1 <- sum(a549[[1]]==0)/(dim(a549[[1]])[1]*dim(a549[[1]])[2])
a2 <- sum(a549[[2]]==0)/(dim(a549[[2]])[1]*dim(a549[[2]])[2])

