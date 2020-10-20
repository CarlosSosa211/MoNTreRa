library(coin) 
set.seed(1)

output <- read.table("../../Carlos/Results/Diff_Ang_432x5Sim_Tissue18/endTreatTumDens.res", header = FALSE)
x <- output[, 1] #no angiogenesis
y <- output[, 2] #angiogenesis

p1 <- hist(x, breaks = 50)                     
p2 <- hist(y, breaks = 50)                     
plot(p1, col = rgb(0, 0, 1, 1/4), xlim = c(0, 30))  
plot(p2, col = rgb(1, 0, 0, 1/4), xlim = c(0, 30), add = T)
p3 <- hist(x - y, breaks = 50)  

wilcox.test(x, y, alternative = "t")
wilcox.test(x, y, paired = TRUE, alternative = "t")
