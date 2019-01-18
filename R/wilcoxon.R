library(coin) 
set.seed(1)

endTreatTumDens <- read.table("../../Carlos/Results/Diff_Ang_432Sim_AllTissues/Tissue4/endTreatTumDens.res", header = FALSE)
x <- endTreatTumDens[, 1] #no angiogenesis
y <- endTreatTumDens[, 2] #angiogenesis

p1 <- hist(x, breaks = 50)                     
p2 <- hist(y, breaks = 50)                     
plot(p1, col = rgb(0, 0, 1, 1/4), xlim = c(0, 30))  
plot(p2, col = rgb(1, 0, 0, 1/4), xlim = c(0, 30), add = T)

wilcox.test(x, y, alternative = "t")
wilcox.test(x, y, paired = TRUE, alternative = "t")
