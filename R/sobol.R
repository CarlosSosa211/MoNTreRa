library(sensitivity)
nFactors <- 5
nRep <- 10000

X1 <- data.frame(matrix(runif(nFactors * nRep), nrow = nRep))
X2 <- data.frame(matrix(runif(nFactors * nRep), nrow = nRep))

Sobol <- sobol2007(model = NULL, X1, X2, nboot = 0)
design <- Sobol$X

write(t(data.matrix(design)), file = "../InputFiles/X.dat",
      ncolumns = nFactors,
      append = FALSE, sep = " ")

write(c(nFactors, nRep), file = "../InputFiles/sobolDim.dat",
      ncolumns = nFactors,
      append = FALSE, sep = " ")

y <- unlist(read.delim("../OutputFiles/Y.res", header = FALSE))
tell(Sobol, y)

print(Sobol)
plot(Sobol)