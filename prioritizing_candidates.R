
#species <- "rat"
#species <- "human"
species <- "mouse"

makePlots <- TRUE
writeTables <- TRUE

if (species == "rat") {
  chrSizesFn <- "~/resources/rat_chromosome_sizes"
  subdir <- "TiPeR_and_SCAP"
  inputStats <- paste(subdir, "/circular-alignments-2021-06-14-", species ,".stats.tsv", sep="")
  filterConfig <- list(mismatchRate=0.007, minCells=1, maxCells=15, minUniqueTerminals=0.3, maxEncompassingRate=0.85, minReads=10)
} else if (species == "human") {
  chrSizesFn <- "~/resources/ucsc/hg38/patch12/chromosomes/chromosome_sizes"
  subdir <- "TiPeR_and_SCAP"
  inputStats <- paste(subdir, "/circular-alignments-2021-06-14-", species ,".stats.tsv", sep="")
  filterConfig <- list(mismatchRate=0.015, minCells=1, maxCells=10, minUniqueTerminals=0.4, maxEncompassingRate=0.80, minReads=10)
} else if (species == "mouse") {
  chrSizesFn <- "~/resources/mm10/chromosomes/chromosome_sizes"
  #subdir <- "mouse-dendrites"
  #inputStats <- paste(subdir, "/circular-alignments-2021-07-22.stats.tsv", sep="")
  #filterConfig <- list(mismatchRate=0.02, minCells=1, maxCells=10, minUniqueTerminals=0.5, maxEncompassingRate=0.80, minReads=8)
  subdir <- "mouse-astrocytes-and-neurons"
  inputStats <- paste(subdir, "/circular-alignments-2022-10-03.stats.tsv", sep="")
  dOutFn <- paste0(subdir, "/stats-and-sel-matrix.tsv")
  filterConfig <- list(mismatchRate=0.02, minCells=2, maxCells=60, minUniqueTerminals=0.5, maxEncompassingRate=0.80, minReads=8)
} else {
  stop("invalid species")
}

d <- read.table(file=inputStats, sep="\t", header=TRUE, as.is=TRUE)

readChrSizes <- function(fn) {
  chrSizes <- read.table(file=fn, sep="\t", header=FALSE, as.is=TRUE)
  names(chrSizes) <- c("chr", "size")
  if (substr(chrSizes$chr[1], 1, 3) != "chr") {
    chrSizes$chr <- paste("chr", chrSizes$chr, sep="")
  }
  chrSizes$linearEnd <- cumsum(chrSizes$size/100)
  chrSizes$linearStart <- chrSizes$linearEnd - chrSizes$size/100
  return(chrSizes)
}

localizePlotFn <- function(prefix) {
  paste(subdir, "/", prefix, "-", species, ".pdf", sep="")
}

chrSizes <- readChrSizes(chrSizesFn)

d$linCoord <- d$refStart/100 + chrSizes$linearStart[match(d$chr, chrSizes$chr)]

convToBin <- function(x) floor(x/10)
d$coordBin <- convToBin(d$linCoord)
d$coordBin <- convToBin(d$linCoord)
chrSizes$binStart <- convToBin(chrSizes$linearStart)
chrSizes$binEnd <- convToBin(chrSizes$linearEnd)
chrSizes$binMiddle <- (chrSizes$binStart + chrSizes$binEnd) / 2

bins <- sort(unique(d$coordBin))

computeHeat <- function(d) {
  sapply(bins, function(i) sum(d$coordBin == i))
}

heat <- computeHeat(d)
sortedHeatHead <- head(sort(heat, decreasing=TRUE), 50)

simHeat <- function(x, n=50) {
  h <- head(sort(x, decreasing=TRUE), n)
  s <- sample(1:n, sum(h), replace=TRUE)
  tbl <- sapply(1:n, function(i) sum(s == i))
  sort(tbl, decreasing=TRUE)
}

simHeatPois <- function(x, n=50) {
  h <- head(sort(x, decreasing=TRUE), n)
  z <- rpois(length(x), lambda=mean(h))
  head(sort(z, decreasing=TRUE), n)
}
  
# This measure of clumpiness suffers from low clumpiness values for very thin data.
clumpiness <- function(x, p=0.8, n=50) {
  h <- head(sort(x, decreasing=TRUE), n)
  w <- p^(0:(n-1))
  h <- h / sum(h)
  sum(h * w)
}

expectedClumpiness <- mean(replicate(100, clumpiness(simHeat(heat))))

# We can get around uneven amounts of data the included and excluded data subsets by
# downsampling the larger set.
compareClumpiness <- function(sel, data, label) {
  d1 <- data[sel,]
  d2 <- data[!sel,]
  if (nrow(d1) < nrow(d2)) {
    d2 <- d2[sample(1:nrow(d2), nrow(d1)),]
  } else {
    d1 <- d1[sample(1:nrow(d1), nrow(d2)),]
  }
  stopifnot(nrow(d1) == nrow(d2))
  c1 <- clumpiness(computeHeat(d1))
  c2 <- clumpiness(computeHeat(d2))
  list(label=label, include=c1, exclude=c2, frac=mean(sel))
}


if (makePlots) {
  pdf(file=localizePlotFn("candidate-density"), height=4, width=20)
  plot(bins, heat, type="l", axes=FALSE)
  axis(2)
  abline(v=chrSizes$binEnd, col="gray")
  par(xpd=NA)
  text(chrSizes$binMiddle, max(heat)*(-0.05), labels=chrSizes$chr, cex=0.7)
  dev.off()
}

if (makePlots) {
  pdf(file=localizePlotFn("mismatch-distribution"), height=5, width=5)
  hist(d$mismatches / d$matches, breaks=50, xlab="mismatch rate", main="", col="gray70")
  dev.off()
}

if (makePlots) {
  pdf(file=localizePlotFn("unique-terminals-fraction-distribution"), height=5, width=5)
  hist(d$uniqueTerminalsFrac, breaks=50, xlab="unique terminals fraction", main="", col="gray70")
  dev.off()
}

if (makePlots) {
  pdf(file=localizePlotFn("cell-count-distribution"), height=5, width=5)
  hist(d$cells, breaks=50, xlab="number of cells", main="", col="gray70")
  dev.off()
}

if (makePlots) {
  pdf(file=localizePlotFn("encompassing-fraction"), height=5, width=5)
  hist(d$encompassing / d$reads, breaks=50, xlab="fraction of reads that are encompassing", main="", col="gray70")
  dev.off()
}

if (makePlots) {
  pdf(file=localizePlotFn("overlapping-fraction"), height=5, width=5)
  hist(d$overlapping / d$aligned, breaks=50, xlab="fraction of aligned bases where mates overlap", main="", col="gray70")
  dev.off()
}

sel1 <- d$excludedReads / d$reads < 1                                # more clumpy
sel2 <- d$mismatches / d$matches < filterConfig$mismatchRate         # less clumpy
sel3 <- filterConfig$minCells < d$cells                              # more clumpy
sel4 <- d$cells <= filterConfig$maxCells                             # less clumpy
sel5 <- filterConfig$minUniqueTerminals < d$uniqueTerminalsFrac      # less clumpy
sel6 <- 0.1 < d$encompassing / d$reads                               # more clumpy
sel7 <- d$encompassing / d$reads < filterConfig$maxEncompassingRate  # less clumpy
sel8 <- filterConfig$minReads <= d$reads                             # more clumpy

if (species == "rat") {
  selOrder <- list("2"=sel2, "5"=sel5, "4"=sel4, "8"=sel8, "7"=sel7, "1"=sel1, "3"=sel3, "6"=sel6)
}
if (species == "human") {
  selOrder <- list("2"=sel2, "4"=sel4, "5"=sel5, "8"=sel8, "7"=sel7, "1"=sel1, "3"=sel3, "6"=sel6)
}
if (species == "mouse") {
  #selOrder <- list("4"=sel4, "3"=sel3, "6"=sel6, "1"=sel1, "2"=sel2, "5"=sel5, "8"=sel8, "7"=sel7)
  selOrder <- list("8"=sel8, "3"=sel3, "4"=sel4, "5"=sel5, "2"=sel2, "1"=sel1, "6"=sel6, "7"=sel7)
}

selMx <- cbind(sel1, sel2, sel3, sel4, sel5, sel6, sel7, sel8)+0
dAndSel <- cbind(d, selMx)
write.table(dAndSel, file=dOutFn, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)


combinedSel <- rep(TRUE, nrow(d))
cOrdering <- c()
nOrdering <- c()
for (i in 1:length(selOrder)) {
  nextSelName <- names(selOrder)[i]
  print(compareClumpiness(selOrder[[nextSelName]][combinedSel], d[combinedSel,], label=names(selOrder)[i]))
  cOrdering <- c(cOrdering, clumpiness(computeHeat(d[combinedSel,])))
  nOrdering <- c(nOrdering, sum(combinedSel))
  combinedSel <- combinedSel & selOrder[[nextSelName]]
}
cOrdering <- c(cOrdering, clumpiness(computeHeat(d[combinedSel,])))
nOrdering <- c(nOrdering, sum(combinedSel))

if (makePlots) {
  pdf(file=localizePlotFn("clumpiness-sequential-filtering"), height=5, width=6.5)
  barplot(cOrdering, xlab="sequential filters applied", ylab="clumpiness",
    names.arg=c("{}", paste("+", names(selOrder), sep="")))
  par(xpd=NA)
  text((0:8)*1.2 + 0.65, cOrdering + 0.015, labels=nOrdering, cex=0.8)
  dev.off()
}

if (makePlots) {
  pdf(file=localizePlotFn("sequential-filtering-counts"), height=5, width=6.5)
  par(xpd=NA)
  barplot(nOrdering, xlab="sequential filters applied", ylab="count",
    names.arg=c("{}", paste("+", names(selOrder), sep="")))
  text((0:8)*1.2 + 0.65, nOrdering + max(nOrdering) * 0.02, labels=nOrdering, cex=0.7)
  dev.off()
}

good <- sel1 & sel2 & sel3 & sel4 & sel5 & sel6 & sel7 & sel8
withoutSel4 <- sel1 & sel2 & sel3 & sel5 & sel6 & sel7 & sel8
dGood <- d[good,]

goodFiles <- paste(paste(dGood$chr, dGood$refStart-1, dGood$refEnd+1, sep="-"), ".aln", sep="")
alnListFn <- paste(subdir, "/good-alignments-", species, ".txt", sep="")

dFiles <- paste(paste(d$chr, d$refStart-1, d$refEnd+1, sep="-"), ".aln", sep="")

if (writeTables) {
  writeLines(goodFiles, con=alnListFn)
}

if (writeTables) {
  goodDataFn <- paste(subdir, "/good-alignments-data-", species, ".tsv", sep="")
  write.table(dGood, file=goodDataFn, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
}

d$overlappingFrac <- d$overlapping / d$aligned
d$mismatchRate <- d$mismatches / (d$matches + d$mismatches)

#expectedHeat <- apply(replicate(100, simHeat(heat)), 1, mean)
#
#
#if (makePlots) {
#  pdf(file=localizePlotFn("hotspot-heat"), height=4.5, width=6)
#  plot(sortedHeatHead, pch=20, xlab="sorted hotspots", ylab="hotspot intensity")
#  lines(1:length(expectedHeat), expectedHeat, col="red")
#  legend("topright", inset=0.03, legend=c("observed hotspot heat", "expected heat under uniform distr."), pch=c(20, NA), col=c("black", "red"), cex=0.9, lty=c(NA,1))
#  dev.off()
#}
#
#coordBinCounts <- table(d$coordBin)
#coordBinCounts <- sort(coordBinCounts, decreasing=TRUE)
#inHotspot <- d$coordBin %in% as.numeric(names(coordBinCounts)[coordBinCounts >= 25])


# END
