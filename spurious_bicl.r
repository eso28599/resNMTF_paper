# removal of spurious biclusters
source("utils.r")
# Libraries
library(philentropy)
library("Matrix")
library("clue")
library("aricode")
library("rio")
library("eList")
library(foreach)
library(doParallel)
library(doSNOW)
library(MASS)
##functions to remove spurious biclusters
###functions for calculating JSD 
jsd_calc <- function(x1, x2){
    #calculate the Jensen Shannon divergence between 
    #vectors x1 and x2
    max_val <- max(x1,x2)
    d1  <- density(x1, from=0, to=max_val)
    d2  <- density(x2, from=0, to=max_val)
    d1$y[d1$x>max(x1)] <- 0
    d2$y[d2$x>max(x2)] <- 0
    return(suppressMessages(JSD(rbind(d1$y, d2$y), unit = "log2", est.prob ="empirical")))
}


check_biclusters <- function(Xinput, Foutput, repeats){
  # calculate JSD between returned F and noise
  # returns scores, avg score and max score
  n_views <- length(Xinput)
  n_clusts <- dim(Foutput[[1]])[2]
  # updated results
  scores <- matrix(0, nrow = n_views, ncol = n_clusts)
  thresholds <- get_thresholds(Xinput, Foutput, repeats)
  for (i in 1:n_views){
    x_noise <- thresholds$data[[i]]
    for (k in 1:n_clusts){
      x <- Foutput[[i]][, k]
      scores[i, k] <- mean(apply(x_noise, 2, 
           function(y) jsd_calc(x,y)))
    }
  }
  return(list("score" = scores,
          "avg_threshold" = thresholds$avg_score,
          "max_threshold" = thresholds$max_score))
}