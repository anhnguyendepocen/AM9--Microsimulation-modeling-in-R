## since we are not doing sex-dependent mortality for microsim Sick-Sicker, we will modify the "mortProb_csv" dataset 
p_mort <- read.csv("mortProb.csv")                     
p_mort <- p_mort[p_mort$Age <= 107,]
mort <- rate <- c()
for (i in seq(1, nrow(p_mort), 2)) {
  mort <- c(mort, (p_mort$p.HD[i] + p_mort$p.HD[i+1])/2)
  rate <- c(rate, (p_mort$r.HD[i] + p_mort$r.HD[i+1])/2)
}
p_mort1 <- data.frame(Age = 0:107, r_HD = rate, p_HD = mort)
write.csv(p_mort1, "mortProb_age.csv", row.names = F)
