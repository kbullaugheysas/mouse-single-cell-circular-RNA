fn <- "/kimdata/bullaugh/analysis/circular-RNA/realignment/mouse-dendrites/breakpoint_donor_A.tsv"
d <- read.table(file=fn, header=TRUE, sep="\t", as.is=TRUE)

chrSizes <- read.table(file="/home/bullaugh/resources/mm10/chromosomes/chromosome_sizes", header=FALSE, sep="\t", as.is=TRUE)
names(chrSizes) <- c("chr", "length")
chrSizes$cummulative <- head(cumsum(c(0, as.numeric(chrSizes$length))), nrow(chrSizes))

# Limit ourselves to the main chromosomes.
d <- d[d$chr %in% chrSizes$chr,]

d$cummulative <- d$pos + chrSizes$cummulative[match(d$chr, chrSizes$chr)]
d$mb <- floor(d$cummulative / 1000000)

h <- data.frame(mb=seq(0, floor(sum(chrSizes$length) / 1000000)))
tbl <- table(d$mb)
h$count <- 0
h$count[match(as.numeric(names(tbl)), h$mb)] <- as.numeric(tbl)

stopifnot(sum(is.na(d$mb)) == 0)
stopifnot(sum(is.na(h$count)) == 0)

pdf(file="mouse-dendrites-hotspots.pdf", height=4, width=15)
plot(h$mb, h$count, type="l", xlab="chromosome", axes=FALSE, ylab="chimeric mappings")
axis(2)
abline(v=c(0, chrSizes$cummulative / 1e6), col="gray")
chrMids <- (chrSizes$cummulative + cumsum(as.numeric(chrSizes$length))) / 2 / 1e6
par(xpd=NA)
text(chrMids, max(h$count)*(-0.07), labels=chrSizes$chr, adj=c(0.5, 0), cex=0.7)
dev.off()

pdf(file="mouse-dendrites-hotspots-log.pdf", height=4, width=15)
plot(h$mb, h$count, type="l", xlab="chromosome", axes=FALSE, ylab="chimeric mappings", log="y")
axis(2)
abline(v=c(0, chrSizes$cummulative / 1e6), col="gray")
chrMids <- (chrSizes$cummulative + cumsum(as.numeric(chrSizes$length))) / 2 / 1e6
par(xpd=NA)
text(chrMids, max(h$count)*(-0.07), labels=chrSizes$chr, adj=c(0.5, 0), cex=0.7)
dev.off()

dChr9 <- d[d$chr == "chr9",]
dChr14 <- d[d$chr == "chr14",]






