doseToControl <- read.table("../../Carlos/Results/TCP/doseToControl_0.res", header = FALSE)
doseToControl <- rnorm(100)
p <- ecdf(doseToControl[,1])
plot(p)
