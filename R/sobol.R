library(sensitivity)
nFactors <- 5
nRep <- 1e5

X1 <- data.frame(matrix(runif(nFactors * nRep), nrow = nRep))
X2 <- data.frame(matrix(runif(nFactors * nRep), nrow = nRep))

Sobol <- sobolmartinez(model = NULL, X1, X2)
design <- Sobol$X

write(t(data.matrix(design)), file = "../InputFiles/X.dat",
      ncolumns = nFactors,
      append = FALSE, sep = " ")

write(c(nFactors, nRep), file = "../InputFiles/sobolDim.dat",
      ncolumns = nFactors,
      append = FALSE, sep = " ")

#--------------------------------------------------------------

yEndTreatTumDens <- unlist(read.delim("../OutputFiles/YEndTreatTumDens.res", header = FALSE))
y3MonTumDens     <- unlist(read.delim("../OutputFiles/Y3MonTumDens.res", header = FALSE))
yRecTumDens      <- unlist(read.delim("../OutputFiles/YRecTumDens.res", header = FALSE))
yTumVol          <- unlist(read.delim("../OutputFiles/YTumVol.res", header = FALSE))
yIntTumDens      <- unlist(read.delim("../OutputFiles/YIntTumDens.res", header = FALSE))
yTimeTo95        <- unlist(read.delim("../OutputFiles/YTimeTo95.res", header = FALSE))
yTimeTo99        <- unlist(read.delim("../OutputFiles/YTimeTo99.res", header = FALSE))
yRecTime         <- unlist(read.delim("../OutputFiles/YRecTime.res", header = FALSE))

tell(Sobol, y)

print(Sobol)
plot(Sobol)