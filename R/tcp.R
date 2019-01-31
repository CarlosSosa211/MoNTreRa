doseToControl2Gy <- read.table("../../Carlos/Results/TCP_0.8_1_0.1_1/doseToControl_2Gy.res", header = FALSE)
doseToControl3Gy <- read.table("../../Carlos/Results/TCP_0.8_1_0.1_1/doseToControl_3Gy.res", header = FALSE)
doseToControl4Gy <- read.table("../../Carlos/Results/TCP_0.8_1_0.1_1/doseToControl_4Gy.res", header = FALSE)

p2 <- ecdf(doseToControl2Gy[, 1])
p3 <- ecdf(doseToControl3Gy[, 1])
p4 <- ecdf(doseToControl4Gy[, 1])

plot(p2, col = rgb(0, 0, 1, 1/4),  xlim = c(0, 100), main = "Tumour control probability", xlab = "Dose (Gy)", ylab = "TCP")
plot(p3, col = rgb(1, 0, 0, 1/4), add = T)
plot(p4, col = rgb(0, 1, 0, 1/4), add = T)
