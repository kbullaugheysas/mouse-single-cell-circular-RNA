configs <- list()
configs[[1]] <- list(prefix="candidates-2020-12-18-coordinates.levenshtein", title="Levenshtein distance between inside and flanking kmers\n(min of both combinations)")
configs[[2]] <- list(prefix="candidates-2020-12-18-coordinates.levenshtein.inside", title="Levenshtein distance between two inside kmers")
configs[[3]] <- list(prefix="candidates-2020-12-18-coordinates.levenshtein.outside", title="Levenshtein distance between two flanking kmers")
configs[[4]] <- list(prefix="candidates-2020-12-18-coordinates.levenshtein.rc", title="Levenshtein distance between inside and flanking kmers (rc)\n(min of both combinations)")
configs[[5]] <- list(prefix="candidates-2020-12-18-coordinates.levenshtein.inside.rc", title="Levenshtein distance between two inside kmers (rc)")
configs[[6]] <- list(prefix="candidates-2020-12-18-coordinates.levenshtein.outside.rc", title="Levenshtein distance between two flanking kmers (rc)")
for (k in c(6, 8, 10, 12, 15, 20, 25)) {
  prefix <- paste("candidates-2020-12-18-coordinates.levenshtein.k", k, sep="")
  title <- paste("Levenshtein distance between inside and flanking kmers (k=", k, ")\n(min of both combinations)", sep="")
  configs[[length(configs)+1]] <- list(prefix=prefix, title=title)
}

plotConfigs <- function(file, configs) {
  pdf(file=file, height=5, width=7)
  for (cfg in configs) {
    if (is.null(cfg$files)) {
      cfg$files <- c(paste(cfg$prefix, ".tsv", sep=""), paste(cfg$prefix, ".null.tsv", sep=""))
    }
    if (is.null(cfg$labels)) {
      cfg$labels <- c("observed distribution", "null distribution")
    }
    d <- read.table(file=cfg$files[1], sep="\t", as.is=TRUE, header=FALSE)
    dNull <- read.table(file=cfg$files[2], sep="\t", as.is=TRUE, header=FALSE)
    d$frac <- d$V2 / sum(d$V2)
    dNull$frac <- dNull$V2 / sum(dNull$V2)
    plot(range(d$V1), range(c(d$frac, dNull$frac)), type="n", xlab="Levenshtein distance", ylab="Probability", main=cfg$title)
    lines(d$V1, d$frac, col="blue")
    points(d$V1, d$frac, pch=21, col="blue", bg="steelblue")
    lines(dNull$V1, dNull$frac, col="gray30")
    points(dNull$V1, dNull$frac, pch=21, col="gray30", bg="gray60")
    legend("topright", legend=cfg$labels, pch=21, col=c("blue", "gray30"), pt.bg=c("steelblue", "gray60"), cex=0.8)
  }
  dev.off()
}

plotConfigs("candidates-2020-12-18-coordinates.levenshtein.pdf", configs)

mouseDendriteConfigs <- list()
mouseDendriteConfigs[[1]] <- list(prefix="mouse-dendrites-candidates-2021-07-22.levenshtein", title="levenshtein distance between inside and flanking kmers\n(min of both combinations); Mouse dendrites; k=10")
mouseDendriteConfigs[[2]] <- list(prefix="mouse-dendrites-candidates-2021-07-22.levenshtein.k8", title="levenshtein distance between inside and flanking kmers\n(min of both combinations); Mouse dendrites; k=8")
mouseDendriteConfigs[[3]] <- list(prefix="mouse-dendrites-candidates-2021-07-22.levenshtein.k12", title="levenshtein distance between inside and flanking kmers\n(min of both combinations); Mouse dendrites; k=12")

plotConfigs("mouse-dendrites-candidates-2021-07-22.levenshtein.pdf", mouseDendriteConfigs)

lowPriorityConfigs <- list()
lowPriorityConfigs[[1]] <- list(prefix="candidates-2020_12_15-low-priority.levenshtein", title="levenshtein distance between inside and flanking kmers\n(min of both combinations); Mouse neurons 2020-12; low priority")

plotConfigs("candidates-2020-12-15-low-priority.pdf", lowPriorityConfigs)

priorityComparison <- list()
priorityComparison[[1]] <- list(files=c("candidates-2020-12-18-coordinates.levenshtein.tsv", "candidates-2020_12_15-low-priority.levenshtein.tsv"),
  title="levenshtein distance between inside and flanking kmers\nMouse neurons 2020-12; comparing priorities", labels=c("high priority", "low priority"))
plotConfigs("candidates-2020-12-15-comparing-priorities.pdf", priorityComparison)

# END
