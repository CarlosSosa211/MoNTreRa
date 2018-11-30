library(sensitivity)
boundaries <- read.delim("../InputFiles/refParIntRT.dat", header = FALSE)
nFactors <- 34
nLevels <- 20
nRep <- 100

Morris <- morris(model = NULL, factors = nFactors, r = nRep,
                 design = list(type = "oat", levels = nLevels, grid.jump = nLevels / 2),
                 binf = boundaries[, 1],
                 bsup = boundaries[, 2],
                 scale = TRUE)
design <- Morris$X

write(t(data.matrix(design)), file = "../InputFiles/X.dat",
      ncolumns = nFactors,
      append = FALSE, sep = " ")

write(c(nFactors, nRep), file = "../InputFiles/morrisDim.dat",
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

tell(Morris, yEndTreatTumDens)

mu.star <- apply(Morris$ee, 2, function(Morris) mean(abs(Morris)))
sigma <- apply(Morris$ee, 2, sd)
print(Morris)
plot(Morris)