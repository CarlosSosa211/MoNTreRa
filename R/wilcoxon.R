library(coin) 
set.seed(1)

endTreatTumDens <- read.table("../../Carlos/Results/Err_Ang_432Sim_AllTissues/Tissue18/endTreatTumDens.res", header = FALSE)
x <- endTreatTumDens[, 1] #no angiogenesis
y <- endTreatTumDens[, 2] #angiogenesis

x <- runif(500, 0.0, 1.0)
y <- x - 0.01 + 0.02 * runif(500, 0.0, 1.0)

x <- rnorm(500, 4.0, 1.0)                     
y <- rnorm(500, 4.0, 1.0)  

x <- c(0.45, 0.5, 0.61, 0.63, 0.75, 0.85, 0.93)
y <- c(0.44, 0.45, 0.52, 0.53, 0.56, 0.58, 0.58, 0.65, 0.79)

p1 <- hist(x, breaks = 50)                     
p2 <- hist(y, breaks = 50)                     
plot(p1, col = rgb(0, 0, 1, 1/4), xlim = c(0, 10))  
plot(p2, col = rgb(1, 0, 0, 1/4), xlim = c(0, 10), add = T)

wilcox.test(x, y, alternative = "t")
wilcox.test(x, y, paired = TRUE, alternative = "t")