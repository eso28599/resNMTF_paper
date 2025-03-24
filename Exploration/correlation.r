#correlation between f score and BiS - phi 4
colnames <- c("Method", "Factor", "Choice", "Measure", "Mean","S.d.")
results_bicl <- read.csv("increasing_bicl/all_results.csv", row.names=1)
results_views <- read.csv("increasing_views2/all_results.csv", row.names=1)
results_indiv <- read.csv("increasing_indiv4/all_results.csv", row.names=1)
results_noise <- read.csv("increasing_noise4/all_results.csv", row.names=1)
names(results_bicl) <- colnames
names(results_views) <- colnames
names(results_indiv) <- colnames
names(results_noise) <- colnames

all_results <- rbind(results_bicl, results_views,
                 results_indiv, results_noise)

f_score <- all_results$Mean[all_results$Measure=="F_score"]
BiS <- all_results$Mean[all_results$Measure=="BiS"]
print(cor(f_score, BiS))


results4 <- read.csv("increasing_phi4/all_results.csv")
results2 <- read.csv("increasing_phi2/all_results.csv")
results3 <- read.csv("increasing_phi3/all_results.csv")
results6 <- read.csv("increasing_phi6/all_results.csv")

all_results <- rbind(results4, results2,
                 results3, results6)
f_score <- all_results$Mean[all_results$Measure=="F_score"]
BiS <- all_results$Mean[all_results$Measure=="BiS"]
print(cor(f_score, BiS))

for(res in list(results4,results2, results3, results6)){
    f_score <- res$Mean[res$Measure=="F_score"]
    BiS <- res$Mean[res$Measure=="BiS"]
    print(cor(f_score, BiS))
}