#correlation between f score and BiS - phi 4
results <- read.csv("sim_folders/increasing_phi4/all_results.csv")
results2 <- read.csv("sim_folders/increasing_phi2/all_results.csv")
results3 <- read.csv("sim_folders/increasing_phi3/all_results.csv")
results6 <- read.csv("sim_folders/increasing_phi6/all_results.csv")
results7 <- read.csv("sim_folders/increasing_phi7/all_results.csv")

all_results <- rbind(results, results2, results3, results6, results7)


f_score <- all_results$Mean[all_results$Measure=="F_score"]
BiS <- all_results$Mean[all_results$Measure=="BiS"]
rel <- all_results$Mean[all_results$Measure=="Rel"]
rec <- all_results$Mean[all_results$Measure=="Rec"]
CSR <- all_results$Mean[all_results$Measure=="CSR"]
Rest <- all_results$Mean[all_results$Measure=="Restrictions"]
cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)
cor(CSR, BiS)
cor(Rest, rec)


#correlation between f score and BiS - phi 3
results <- read.csv("increasing_phi3/all_results.csv")
f_score <- results$Mean[results$Measure=="F_score"]
BiS <- results$Mean[results$Measure=="BiS"]
rel <- results$Mean[results$Measure=="Rel"]
rec <- results$Mean[results$Measure=="Rec"]
cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)


#correlation between f score and BiS - phi 3
results <- read.csv("increasing_phi3/all_results.csv")
f_score <- results$Mean[results$Measure=="F_score"]
BiS <- results$Mean[results$Measure=="BiS"]
rel <- results$Mean[results$Measure=="Rel"]
rec <- results$Mean[results$Measure=="Rec"]
cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)
cor(rec, rel)

#correlation between f score and BiS - phi 2
results <- read.csv("increasing_phi2/all_results.csv")
f_score <- results$Mean[results$Measure=="F_score"]
BiS <- results$Mean[results$Measure=="BiS"]
rel <- results$Mean[results$Measure=="Rel"]
rec <- results$Mean[results$Measure=="Rec"]
cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)
cor(rec, rel)


#correlation between f score and new BiS
results <- read.csv("results/all_results.csv")
f_score <- results$Mean[results$Measure=="F_score"]
BiS <- results$Mean[results$Measure=="BiS"]
rel <- results$Mean[results$Measure=="Rel"]
rec <- results$Mean[results$Measure=="Rec"]
cor(f_score, BiS)
cor(rel, BiS)
cor(rec, BiS)
